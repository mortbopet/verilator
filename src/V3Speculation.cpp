// -*- mode: C++; c-file-style: "cc-mode" -*-
//*************************************************************************
// DESCRIPTION: Verilator: Constant folding
//
// Code available from: https://verilator.org
//
//*************************************************************************
//
// Copyright 2003-2021 by Wilson Snyder. This program is free software; you
// can redistribute it and/or modify it under the terms of either the GNU
// Lesser General Public License Version 3 or the Perl Artistic License
// Version 2.0.
// SPDX-License-Identifier: LGPL-3.0-only OR Artistic-2.0
//
//*************************************************************************

#include "V3Speculation.h"
#include "V3Ast.h"
#include "V3PartitionGraph.h"
#include "V3AstDot.h"
#include "V3Dataflow.h"

#include <algorithm>

//######################################################################
// SpeculationReplaceVisitor

inline string specNameSuffix(bool branch) { return "__SPEC__" + string(branch ? "t" : "f"); }

class SpeculativeReplaceVisitor final : public AstNVisitor {
private:
    // NODE STATE
    // AstNode::user()  -> node was visited
    AstUser1InUse m_inuser1;

    // Mapping between non-speculative to replaced speculative variable, which is being written to
    // in the partition.
    VarReplMapping m_specOutVars;
    AstVar* m_varp;  // Boolean variable to speculate on
    bool m_branch;
    AstMTaskBody* m_mtaskp;
    AstNodeModule* m_modp;
    int it = 0;

    string nameSuffix() const { return specNameSuffix(m_branch); }

    void dumpAST() {
        AstDotDumper dump(m_mtaskp);
        dump.dumpDotFilePrefixedAlways(cvtToStr(m_mtaskp->name() + "_" + cvtToStr(it++)));
    }

    void replaceBool(AstNode* oldp, bool value) {
        UASSERT(oldp, "Null old");
        AstNode* newp = new AstConst(oldp->fileline(), value);
        newp->dtypeFrom(oldp);
        oldp->replaceWith(newp);
        VL_DO_DANGLING(oldp->deleteTree(), oldp);
    }

    /**
     * @brief findModp
     * recursively backtrack the parent of @p nodep to locate the parent module.
     */
    AstModule* findModp(AstNode* nodep) {
        while (nodep != nullptr) {
            auto* modp = dynamic_cast<AstModule*>(nodep);
            if (modp) {
                return modp;
            } else {
                nodep = nodep->backp();
            }
        }
        UASSERT(nodep, "Tried to find parent module of node not inside a module");
        return nullptr;
    }

    void visit(AstNode* nodep) {
        assert(nodep->user1() == 0 && "!");

        nodep->user1(1);
        dumpAST();
        iterateChildren(nodep);
    }

    void visit(AstNodeAssign* assignp) {
        auto* varrefp = dynamic_cast<AstVarRef*>(assignp->lhsp());
        UASSERT(varrefp, "Assignment to non-variable?");
        AstVar* prevarp = varrefp->varp();
        UASSERT(prevarp != m_varp, "Assignment to variable we speculate on?");

        assignp->dump(cout);
        cout << endl;

        // Find (or construct) speculative output variable
        auto replVarp = m_specOutVars.find(prevarp);
        if (replVarp == m_specOutVars.end()) {
            AstVar* specVarp = prevarp->cloneTree(false);
            specVarp->name(specVarp->name() + nameSuffix());
            m_modp->addStmtp(specVarp);
            replVarp = m_specOutVars.insert(m_specOutVars.end(),
                                            {prevarp, VarReplInfo{assignp, varrefp, specVarp}});
        }

        // And replace
        varrefp->varp(replVarp->second.replvarp);

        visit(static_cast<AstNode*>(assignp));
    }

    void visit(AstVarRef* varrefp) {
        AstVar* varp = varrefp->varp();
        if (varp == m_varp) {
            replaceBool(varrefp, m_branch);
            return;
        }

        auto it = m_specOutVars.find(varp);
        if (it != m_specOutVars.end()) {
            // Replace reference with speculative variable
            varrefp->varp(it->second.replvarp);
        }

        visit(static_cast<AstNode*>(varrefp));
    }

    void visit(AstCCall* ccallp) {
        if (ccallp->user1()) return;
        ccallp->user1(1);

        // Create a duplicate of the called C function, and perform boolean const propagation
        // into this function.
        AstCFunc* specfuncp = ccallp->funcp()->cloneTree(false);
        specfuncp->name(specfuncp->name() + nameSuffix());
        m_modp->addStmtp(specfuncp);

        // Modify target function call to the new, speculative function
        ccallp->funcp(specfuncp);

        // Iterate into the speculative function, further replacing occurances of m_varp
        visit(specfuncp);
    }

public:
    explicit SpeculativeReplaceVisitor(AstNodeModule* modp, AstMTaskBody* bodyp, AstVar* varp,
                                       bool branch)
        : m_varp(varp)
        , m_branch(branch)
        , m_mtaskp(bodyp)
        , m_modp(modp) {
        visit(bodyp);
    }

    const VarReplMapping& specOutVars() const { return m_specOutVars; }
};

//######################################################################
// V3Speculation

V3Speculation::V3Speculation() {
    // Process each module in turn
    for (AstNodeModule* nodep = v3Global.rootp()->modulesp(); nodep;
         nodep = VN_CAST(nodep->nextp(), NodeModule)) {
        speculateModule(nodep);
    }
}

class GatherScopeNamesVisitor final : public AstNVisitor,
                                      public std::unordered_map<AstVar*, std::string> {
private:
    void visit(AstNode* nodep) { iterateChildren(nodep); }
    void visit(AstNodeVarRef* varrefp) {
        (*this)[varrefp->varp()] = varrefp->hiernameToProt();
        visit(static_cast<AstNode*>(varrefp));
    }
    void visit(AstCCall* ccallp) { visit(ccallp->funcp()); }

public:
    GatherScopeNamesVisitor(ExecMTask* mtaskp) { iterate(mtaskp->bodyp()); }
};

void V3Speculation::doSpeculation(AstNodeModule* modp, const Speculateable& s) {
    auto* execGraphp = v3Global.rootp()->execGraphp();
    updateDataflowInfo(modp);

    GatherScopeNamesVisitor scopeNames(s.cons);

    m_dfgs[s.cons]->dumpDotFilePrefixedAlways("DFG_" + cvtToStr(s.cons->id()));
    AstDotDumper dump(s.cons->bodyp());
    dump.dumpDotFilePrefixedAlways(cvtToStr(s.cons->id()));

    // Duplicate the MTask, representing the true and false constant value of the speculation
    // variable
    auto* consmtbodyp_t = s.cons->bodyp()->cloneTree(false);
    auto* consmtbodyp_f = s.cons->bodyp()->cloneTree(false);

    AstDotDumper dump3(consmtbodyp_t);
    dump3.dumpDotFilePrefixedAlways("Truebody_copy_pre");

    // Speculate into them. Manage scope of SpeculativeReplaceVisitors to ensure destruction of
    // user pointers.
    auto performSpeculation = [&](AstMTaskBody* mtaskbody, bool branch) {
        SpeculativeReplaceVisitor specVisitor(modp, mtaskbody, s.specVar, branch);
        return specVisitor.specOutVars();
    };
    VarReplMapping specTOutVars = performSpeculation(consmtbodyp_t, true);
    VarReplMapping specFOutVars = performSpeculation(consmtbodyp_f, false);

    // Insert them into the dependency graph
    execGraphp->addMTaskBody(consmtbodyp_t);
    execGraphp->addMTaskBody(consmtbodyp_f);

    // Create ExecMTasks
    auto* consmtp_t
        = new ExecMTask(execGraphp->mutableDepGraphp(), consmtbodyp_t, m_nextMTaskID++);
    auto* consmtp_f
        = new ExecMTask(execGraphp->mutableDepGraphp(), consmtbodyp_f, m_nextMTaskID++);
    consmtp_t->speculative(ExecMTask::Speculative::True, s.cons, consmtp_f);
    consmtp_f->speculative(ExecMTask::Speculative::False, s.cons, consmtp_t);
    consmtbodyp_t->execMTaskp(consmtp_t);
    consmtbodyp_f->execMTaskp(consmtp_f);

    // Create dependency edges for all incoming variables except the speculated boolean. We cannot
    // use the dependency graph here, due to transitive edge removal, and by the fact that the only
    // remaining dependency edge in the graph is the actual boolean which we are currently
    // speculating on.
    // @todo: Only create edge if actual dependency (ie. avoid loops, we need topological sort)
    for (const auto& invarp : m_io.at(modp).at(s.cons).ins) {
        if (invarp == s.specVar) continue;
        if (invarp->varType() == AstVarType::PORT) continue;

        ExecMTask* prodMTaskp = m_varProducedBy.at(invarp);
        if (prodMTaskp == s.cons) { continue; }
        assert(prodMTaskp != consmtp_t);
        new V3GraphEdge(execGraphp->mutableDepGraphp(), prodMTaskp, consmtp_t, 1);
        assert(prodMTaskp != consmtp_f);
        new V3GraphEdge(execGraphp->mutableDepGraphp(), prodMTaskp, consmtp_f, 1);
    }

    // Add speculative resolution statements to end of speculative bodies
    auto* condp_t = new AstVarRef(s.specVar->fileline(), s.specVar, VAccess::READ);
    condp_t->hiernameToProt(scopeNames.at(s.specVar));
    consmtbodyp_t->addStmtsp(new AstSpecResolveBool(s.cons->bodyp()->fileline(), condp_t,
                                                    genCommitSpecVarStmts(specTOutVars), s.prod));

    auto* condp_f = new AstVarRef(s.specVar->fileline(), s.specVar, VAccess::READ);
    condp_f->hiernameToProt(scopeNames.at(s.specVar));
    consmtbodyp_f->addStmtsp(new AstSpecResolveBool(s.cons->bodyp()->fileline(),
                                                    new AstNot(condp_f->fileline(), condp_f),
                                                    genCommitSpecVarStmts(specFOutVars), s.prod));

    // The speculated mtasks should inherit the outgoing edges from their origin MTask
    for (V3GraphEdge* edgep = s.cons->outBeginp(); edgep; edgep = edgep->outNextp()) {
        assert(consmtp_t != edgep->top());
        new V3GraphEdge(execGraphp->mutableDepGraphp(), consmtp_t, edgep->top(), 1);
        assert(consmtp_f != edgep->top());
        new V3GraphEdge(execGraphp->mutableDepGraphp(), consmtp_f, edgep->top(), 1);
    }

    // Finally, delete the original MTask
    s.cons->bodyp()->unlinkFrBack()->deleteTree(), s.cons->bodyp();
    s.cons->unlinkDelete(execGraphp->mutableDepGraphp());

    AstDotDumper dump2(modp);
    dump2.dumpDotFilePrefixedAlways(modp->name() + "ast_post_spec");
    execGraphp->mutableDepGraphp()->dumpDotFilePrefixedAlways("dep_post_spec");
}

AstNode* V3Speculation::genCommitSpecVarStmts(const VarReplMapping& mapping) {
    AstNode* stmtsp = nullptr;
    for (const auto& it : mapping) {
        // Base varrefs on original var reference (maintains scope, fileline info etc.)
        auto* fromp = it.second.varrefp->cloneTree(false);
        auto* top = it.second.varrefp->cloneTree(false);

        fromp->varp(it.second.replvarp);
        fromp->access(VAccess::READ);
        top->varp(it.first);
        fromp->access(VAccess::WRITE);

        auto* assignp = new AstAssign(it.second.assignp->fileline(), top, fromp);

        stmtsp = AstNode::addNext(stmtsp, assignp);
    }
    return stmtsp;
}

void V3Speculation::speculateModule(AstNodeModule* modp) {
    V3Graph* mtasksp = v3Global.rootp()->execGraphp()->mutableDepGraphp();

    // Locate module variables
    std::vector<AstVar*> boolVars1;
    for (AstNode* nodep = modp->stmtsp(); nodep; nodep = nodep->nextp()) {
        if (AstVar* varp = dynamic_cast<AstVar*>(nodep)) {
            // Reset producer/consumer IDs (updated later)
            varp->clearMTaskIds();
            if (varp->widthMin() == 1) { boolVars1.push_back(varp); }
        }
    }

    updateDataflowInfo(modp);

    // Dependency filtering
    // Filter to the set of dependencies which result as an in-evaluation dependency
    // rank(producing node) < rank(consuming node)
    mtasksp->order();
    std::vector<AstVar*> boolVars2;
    std::copy_if(boolVars1.begin(), boolVars1.end(), std::back_inserter(boolVars2),
                 [this](const AstVar* varp) {
                     if (varp->prodMtaskIds().size() == 0) { return false; }
                     const auto& prodBy = m_mtaskIdToMTask.find(*varp->prodMtaskIds().begin());
                     if (prodBy == m_mtaskIdToMTask.end()) {
                         // Produced by other MTask not in this module, skip...
                         return false;
                     }

                     // Produced by
                     const int prodRank = prodBy->second->rank();

                     return std::any_of(
                         varp->consMtaskIds().begin(), varp->consMtaskIds().end(),
                         [&](int mtaskid) {
                             const int consRank = m_mtaskIdToMTask.at(mtaskid)->rank();
                             return prodRank < consRank;  // todo: <= here indicates mtasks
                                                          // that might be split-able!
                         });
                 });

    // Locate boolean edges in the dependency graph
    std::map<ExecMTask*, std::vector<Speculateable>> speculateable;

    for (AstVar* varp : boolVars2) {
        ExecMTask* prodp = m_mtaskIdToMTask.at(*varp->prodMtaskIds().begin());

        // Locate MTasks which depend on this variable:
        // - There is an edge between the two MTasks
        // - The given variable is the only depending variable produced by the source MTask
        //
        // Due to transitive edges, the consuming MTask can inherit the incoming dependencies of
        // the producing mtask through speculation.
        for (V3GraphEdge* edgep = prodp->outBeginp(); edgep; edgep = edgep->outNextp()) {
            ExecMTask* consp = dynamic_cast<ExecMTask*>(edgep->top());

            std::vector<AstVar*> sharedVars;

            auto& prodOuts = m_dfgs.at(prodp)->outs();
            auto& consIns = m_dfgs.at(consp)->ins();

            assert(std::find(prodOuts.begin(), prodOuts.end(), varp) != prodOuts.end()
                   && "varp should be in produced set of producing mtask");

            std::set_intersection(prodOuts.begin(), prodOuts.end(), consIns.begin(), consIns.end(),
                                  std::back_inserter(sharedVars));

            if (sharedVars.size() == 1
                && std::find(sharedVars.begin(), sharedVars.end(), varp) != sharedVars.end()) {
                UASSERT(std::find(consIns.begin(), consIns.end(), varp) != consIns.end(),
                        "varp should be in consumed set of consuming mtask");

                Speculateable s;
                s.prod = prodp;
                s.cons = consp;
                s.specVar = varp;
                speculateable[s.cons].push_back(s);
            }
        }
    }

    // For now, only speculate a single variable in partitions which could in theory speculate
    // multiple
    for (auto& it : speculateable) {
        if (it.second.size() > 1) { it.second = {it.second.at(0)}; }
    }

    // Go speculate!
    for (const auto& s : speculateable) { doSpeculation(modp, s.second.at(0)); }

    // New mtasks were inserted, reassure order of graph
    mtasksp->order();
}

void V3Speculation::updateDataflowInfo(AstNodeModule* modp) {
    m_dfgs.clear();
    m_io.clear();
    m_mtaskIdToMTask.clear();
    m_varProducedBy.clear();

    std::vector<AstVar*> variables;
    for (AstNode* nodep = modp->stmtsp(); nodep; nodep = nodep->nextp()) {
        if (AstVar* varp = dynamic_cast<AstVar*>(nodep)) {
            // Reset producer/consumer IDs (updated later)
            varp->clearMTaskIds();
            variables.push_back(varp);
        }
    }

    // Gather DFGs for each MTask. This also updates the producer/consumer info of each variable
    V3Graph* mtasksp = v3Global.rootp()->execGraphp()->mutableDepGraphp();
    for (V3GraphVertex* vxp = mtasksp->verticesBeginp(); vxp; vxp = vxp->verticesNextp()) {
        ExecMTask* mtaskp = dynamic_cast<ExecMTask*>(vxp);
        m_nextMTaskID = m_nextMTaskID <= mtaskp->id() ? mtaskp->id() + 1 : m_nextMTaskID;
        m_dfgs[mtaskp] = new DFG(mtaskp->bodyp());
        m_mtaskIdToMTask[mtaskp->id()] = mtaskp;
    }

    for (const auto& varp : variables) {
        for (const int consMTaskid : varp->consMtaskIds()) {
            auto it = m_mtaskIdToMTask.find(consMTaskid);
            assert(it != m_mtaskIdToMTask.end());
            m_io[modp][it->second].ins.insert(varp);
        }
        for (const int prodMTaskid : varp->prodMtaskIds()) {
            auto it = m_mtaskIdToMTask.find(prodMTaskid);
            assert(it != m_mtaskIdToMTask.end());
            ExecMTask* prodp = it->second;
            m_io[modp][prodp].outs.insert(varp);

            auto preProdBy = m_varProducedBy.find(varp);
            if (preProdBy != m_varProducedBy.end()) {
                bool ok = false;
                // Only valid if speculated
                if (preProdBy->second->speculative() != ExecMTask::Speculative::None
                    && prodp->speculative() != ExecMTask::Speculative::None) {
                    if (preProdBy->second->partnerSpecMTaskp() != prodp->partnerSpecMTaskp()) {
                        ok = true;
                    }
                }
                assert(ok && "Multiple producers only allowed for speculative pairs");
            }
            auto* mtask = m_mtaskIdToMTask.at(prodMTaskid);
            assert(mtask->id() == prodMTaskid
                   && "Mismatched mtask ID and actual produced mtask ID");
            m_varProducedBy[varp] = mtask;
        }
    }

    cout << "IO is:" << endl;
    for (const auto& it : m_io.at(modp)) {
        cout << it.first->name();
        cout << "\t ins: ";
        for (const auto& it2 : it.second.ins) { cout << it2->name(); }
        cout << endl << "\t outs: ";
        for (const auto& it2 : it.second.outs) { cout << it2->name(); }
        cout << endl;
    }
}
