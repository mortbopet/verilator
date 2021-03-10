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
            replVarp = m_specOutVars.insert(m_specOutVars.end(), {prevarp, specVarp});
        }

        // And replace
        varrefp->varp(replVarp->second);

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
            varrefp->varp(it->second);
        }

        visit(static_cast<AstNode*>(varrefp));
    }

    void visit(AstCCall* ccallp) {
        if (ccallp->user1()) return;
        ccallp->user1(1);

        // Create a duplicate of the called C function, and perform boolean const propagation
        // into this function.
        const AstCFunc* origfuncp = ccallp->funcp();
        AstCFunc* specfuncp = ccallp->funcp()->cloneTree(false);  // i think true here?
        specfuncp->name(specfuncp->name() + nameSuffix());
        origfuncp->scopep()->addActivep(specfuncp);

        // Replace the function call to the new, speculative function
        auto newFuncCall = new AstCCall(ccallp->fileline(), specfuncp, ccallp->argsp());
        ccallp->replaceWith(newFuncCall);
        newFuncCall->user1(1);  // prevent visitor recursion into new function call

        dumpAST();

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

void V3Speculation::doSpeculation(AstNodeModule* modp, const Speculateable& s) {
    auto* execGraphp = v3Global.rootp()->execGraphp();

    for (const auto& consmtaskp : s.cons) {
        m_dfgs[consmtaskp]->dumpDotFilePrefixedAlways("DFG_" + cvtToStr(consmtaskp->id()));
        AstDotDumper dump(consmtaskp->bodyp());
        dump.dumpDotFilePrefixedAlways(cvtToStr(consmtaskp->id()));

        // Duplicate the MTask, representing the true and false constant value of the speculation
        // variable
        auto* consmtbodyp_t = consmtaskp->bodyp()->cloneTree(false);
        auto* consmtbodyp_f = consmtaskp->bodyp()->cloneTree(false);

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
        consmtp_t->speculative(ExecMTask::Speculative::True, consmtaskp);
        auto* consmtp_f
            = new ExecMTask(execGraphp->mutableDepGraphp(), consmtbodyp_f, m_nextMTaskID++);
        consmtp_f->speculative(ExecMTask::Speculative::False, consmtaskp);

        // Create dependency edges for all incoming variables except the speculated boolean
        // @todo: Only create edge if actual dependency (ie. avoid loops, we need topological sort)
        for (const auto& in : m_io.at(modp).at(consmtaskp).ins) {
            if (in == s.specVar) continue;
            auto& prodmtasks = in->prodMtaskIds();
            if (prodmtasks.size() == 0) { continue; }
            for (int prodMTaskId : prodmtasks) {
                ExecMTask* prodMTaskp = m_mtaskIdToMTask[prodMTaskId];
                if (prodMTaskp == consmtaskp) { continue; }
                new V3GraphEdge(execGraphp->mutableDepGraphp(), prodMTaskp, consmtp_t, 1);
                new V3GraphEdge(execGraphp->mutableDepGraphp(), prodMTaskp, consmtp_f, 1);
            }
        }

        // Replace original MTask with speculative resolve node
        auto* specResolvep = new AstSpecResolveBool(
            consmtaskp->bodyp()->fileline(),
            new AstVarRef(consmtaskp->bodyp()->fileline(), s.specVar, VAccess::READ), nullptr,
            nullptr);
        consmtaskp->bodyp()->stmtsp()->replaceWith(specResolvep);

        // Create true and false branches for committing the speculative calculated out values
        genCommitSpecVarStmts(specResolvep->ifsp(), specTOutVars);
        genCommitSpecVarStmts(specResolvep->elsesp(), specFOutVars);
        specResolvep->mtaskDepTrue(consmtp_t);
        specResolvep->mtaskDepFalse(consmtp_f);
    }

    AstDotDumper dump2(modp);
    dump2.dumpDotFilePrefixedAlways(modp->name() + "_post_spec");
}

void V3Speculation::genCommitSpecVarStmts(AstNode* stmtsp, const VarReplMapping& mapping) {
    for (const auto& it : mapping) {
        AstNode* fromp = new AstVarRef(it.first->fileline(), it.second, VAccess::READ);
        AstNode* top = new AstVarRef(it.first->fileline(), it.first, VAccess::WRITE);
        AstNode* assignp = new AstAssign(it.first->fileline(), top, fromp);
        stmtsp = AstNode::addNext(stmtsp, assignp);
    }
}

void V3Speculation::speculateModule(AstNodeModule* modp) {
    // Locate module variables
    std::vector<AstVar*> variables;
    for (AstNode* nodep = modp->stmtsp(); nodep; nodep = nodep->nextp()) {
        if (AstVar* varp = dynamic_cast<AstVar*>(nodep)) {
            // Reset producer/consumer IDs (updated later)
            // varp->clearMTaskIds();
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

        AstDotDumper dump(mtaskp->bodyp());
        dump.dumpDotFilePrefixedAlways(cvtToStr(mtaskp->id()));
    }

    gatherIO(modp);

    // Filter the set of boolean variables which are available for speculation.
    // This helps us narrow the search for potential MTasks to speculate.

    // Local filtering
    std::vector<AstVar*> boolVars1;
    std::copy_if(
        variables.begin(), variables.end(), std::back_inserter(boolVars1), [](const AstVar* varp) {
            bool relevant = true;
            relevant &= varp->widthMin() == 1;
            relevant &= varp->consMtaskIds().size() > 0;  // Must be a dependency
            relevant &= varp->prodMtaskIds().size() == 1;  // Must be an output of another Mtask
            return relevant;
        });

    // Dependency filtering
    // Filter to the set of dependencies which result as an in-evaluation dependency
    // rank(producing node) < rank(consuming node)
    mtasksp->order();
    std::vector<AstVar*> boolVars2;
    std::copy_if(
        boolVars1.begin(), boolVars1.end(), std::back_inserter(boolVars2),
        [this](const AstVar* varp) {
            // Produced by
            const int prodRank = m_mtaskIdToMTask.at(*varp->prodMtaskIds().begin())->rank();

            return std::any_of(varp->consMtaskIds().begin(), varp->consMtaskIds().end(),
                               [&](int mtaskid) {
                                   const int consRank = m_mtaskIdToMTask.at(mtaskid)->rank();
                                   return prodRank < consRank;  // todo: <= here indicates mtasks
                                                                // that might be split-able!
                               });
        });

    // Locate boolean edges in the dependency graph
    std::map<AstVar*, Speculateable> speculateable;
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

            UASSERT(std::find(prodOuts.begin(), prodOuts.end(), varp) != prodOuts.end(),
                    "varp should be in produced set of producing mtask");

            std::set_intersection(prodOuts.begin(), prodOuts.end(), consIns.begin(), consIns.end(),
                                  std::back_inserter(sharedVars));

            if (sharedVars.size() == 1
                && std::find(sharedVars.begin(), sharedVars.end(), varp) != sharedVars.end()) {
                UASSERT(std::find(consIns.begin(), consIns.end(), varp) != consIns.end(),
                        "varp should be in consumed set of consuming mtask");

                auto it = speculateable.find(varp);
                if (it == speculateable.end()) {
                    it = speculateable.insert(it, {varp, Speculateable()});
                    it->second.prod = prodp;
                    it->second.specVar = varp;
                }
                UASSERT(it->second.prod == nullptr || it->second.prod == prodp,
                        "Multiple producers not allowed at this point");
                it->second.cons.push_back(consp);
            }
        }
    }

    if (debug() || true) {
        std::cout << "Speculateable variables/partitions are" << endl;
        for (const auto& it : speculateable) {
            std::cout << "\t" << it.first->origName() << " (mt" << it.second.prod->id() << " "
                      << it.second.prod->cost() << ") -> ";
            for (const auto& consp : it.second.cons) {
                std::cout << "(mt" << consp->id() << " " << consp->cost() << ") ";
            }
            std::cout << endl;
        }
    }

    // Go speculate!
    for (const auto& s : speculateable) { doSpeculation(modp, s.second); }

    // New mtasks were inserted, reassure order of graph
    mtasksp->order();
}

void V3Speculation::gatherIO(AstNodeModule* modp) {
    std::vector<AstVar*> variables;
    for (AstNode* nodep = modp->stmtsp(); nodep; nodep = nodep->nextp()) {
        if (AstVar* varp = dynamic_cast<AstVar*>(nodep)) {
            // Reset producer/consumer IDs (updated later)
            // varp->clearMTaskIds();
            variables.push_back(varp);
        }
    }

    for (const auto& varp : variables) {
        for (const int consMTaskid : varp->consMtaskIds()) {
            auto it = m_mtaskIdToMTask.find(consMTaskid);
            if (it == m_mtaskIdToMTask.end()) continue;  // stale
            m_io[modp][it->second].ins.insert(varp);
        }
        for (const int prodMTaskid : varp->prodMtaskIds()) {
            auto it = m_mtaskIdToMTask.find(prodMTaskid);
            if (it == m_mtaskIdToMTask.end()) continue;  // stale
            m_io[modp][it->second].outs.insert(varp);
        }
    }
}
