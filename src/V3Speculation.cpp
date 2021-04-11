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
#include "V3Const.h"

#include <algorithm>

//######################################################################
// SpeculationReplaceVisitor

inline string specNameSuffix(bool branch) { return "__SPEC__" + string(branch ? "t" : "f"); }

inline void replaceBool(AstNode* oldp, bool value) {
    UASSERT(oldp, "Null old");
    AstNode* newp = new AstConst(oldp->fileline(), value);
    newp->dtypeFrom(oldp);
    oldp->replaceWith(newp);
    VL_DO_DANGLING(oldp->deleteTree(), oldp);
}

inline AstVar* replaceVar(AstNode* oldp, string prefix) {
    UASSERT(oldp, "Null old");
    AstVar* dummyvar
        = new AstVar(oldp->fileline(), AstVarType::BLOCKTEMP, prefix, VFlagLogicPacked(), 0);
    AstNode* newp = new AstVarRef(oldp->fileline(), dummyvar, VAccess::READ);
    newp->dtypeFrom(oldp);
    oldp->replaceWith(newp);
    VL_DO_DANGLING(oldp->deleteTree(), oldp);
    return dummyvar;
}

class SpecCondExprReplaceVisitor final : public AstNVisitor {
private:
    const Speculateable& m_s;  // variables/expressions to replace
    AstVar* m_dummyVar = nullptr;
    AstNode* m_oldCondExpr = nullptr;

public:
    void visit(AstNode* nodep) { iterateChildren(nodep); }

    void checkCondp(AstNode* nodep, AstNode* condp) {
        if (nodep == m_s.condSpec.brExpr) {
            // First we clone the conditional expression whereafter it is replaced with a reference
            // to a dummy variable. This is quite redundant, but it gives us a way to track the
            // location of the selected conditional expression, after we clone the tree. We can't
            // use user pointers (not cloned, not allowed by Verilator).
            m_oldCondExpr = condp->cloneTree(false);
            m_dummyVar = replaceVar(condp, "__SPEC__cond_dummy");
        } else {
            visit(static_cast<AstNode*>(nodep));
        }
    }

    void visit(AstNodeCond* nodep) { checkCondp(nodep, nodep->condp()); }
    void visit(AstNodeIf* nodep) { checkCondp(nodep, nodep->condp()); }
    void visit(AstCCall* ccallp) { visit(ccallp->funcp()); }

public:
    explicit SpecCondExprReplaceVisitor(AstMTaskBody* bodyp, const Speculateable& s)
        : m_s(s) {
        visit(bodyp);
    }

    AstVar* dummyVar() const { return m_dummyVar; }
    AstNode* condExpr() const { return m_oldCondExpr; }
};

class SpeculativeReplaceVisitor final : public AstNVisitor {
private:
    // NODE STATE
    // AstNode::user()  -> node was visited
    AstUser1InUse m_inuser1;

    // Mapping between non-speculative to replaced speculative variable, which is being written to
    // in the partition.
    VarReplMapping m_specOutVars;
    Speculateable& m_s;  // variables/expressions to replace
    bool m_branch;
    AstMTaskBody* m_mtaskp;
    AstNodeModule* m_modp;
    int it = 0;

    // Maintain set of replaced functions. When deleting the speculated mtasks' tree, we are not
    // deleting the functions that we've replaced as well, so we have to track them separately.
    std::set<AstCFunc*> m_replacedFunctions;

    string nameSuffix() const { return specNameSuffix(m_branch); }

    void dumpAST() {
        AstDotDumper dump(m_mtaskp);
        dump.dumpDotFilePrefixedAlways(cvtToStr(m_mtaskp->name() + "_" + cvtToStr(it++)));
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

        iterateChildren(nodep);
    }

    void visit(AstNodeAssign* assignp) {
        auto* varrefp = dynamic_cast<AstVarRef*>(assignp->lhsp());
        UASSERT(varrefp, "Assignment to non-variable?");
        AstVar* prevarp = varrefp->varp();
        if (m_s.specBoolVar) {
            UASSERT(prevarp != m_s.specBoolVar, "Assignment to variable we speculate on?");
        }

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

    bool isSpeculatedCondExpr(AstNode* nodep) const {
        if (m_s.condSpec.dummyVar == nullptr) { return false; }
        if (auto* varrefp = dynamic_cast<AstVarRef*>(nodep)) {
            return varrefp->varp() == m_s.condSpec.dummyVar;
        }
        return false;
    }

    // AstNodeCond and AstNodeIf don't inherit from the same parent, but have identical structure,
    // hmpfh...
    void handleBranch(AstNode* parent, AstNode* condp, AstNode* thenp, AstNode* elsep) {
        if (isSpeculatedCondExpr(condp)) {
            // Exchange parent with the branch corresponding to m_branch
            parent->replaceWith(m_branch ? thenp->cloneTree(false) : elsep->cloneTree(false));
            VL_DO_DANGLING(parent->deleteTree(), parent);
        } else {
            if (condp) iterate(condp);
            if (thenp) iterate(thenp);
            if (elsep) iterate(elsep);
        }
    }

    void visit(AstNodeCond* nodep) {
        handleBranch(nodep, nodep->condp(), nodep->expr1p(), nodep->expr2p());
    }

    void visit(AstNodeIf* nodep) {
        handleBranch(nodep, nodep->condp(), nodep->ifsp(), nodep->elsesp());
    }

    void visit(AstNodeVarRef* varrefp) {
        AstVar* varp = varrefp->varp();

        // Only replace when we are not replacing a conditional expression
        if (m_s.specBoolVar && varp == m_s.specBoolVar) {
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

    void visit(AstNodeCCall* ccallp) {
        if (ccallp->user1()) return;
        ccallp->user1(1);

        // Create a duplicate of the called C function, and perform boolean const propagation
        // into this function.
        m_replacedFunctions.insert(ccallp->funcp());
        AstCFunc* specfuncp = ccallp->funcp()->cloneTree(false);
        specfuncp->name(specfuncp->name() + nameSuffix());
        m_modp->addStmtp(specfuncp);

        // Modify target function call to the new, speculative function
        ccallp->funcp(specfuncp);

        // Iterate into the speculative function, further replacing occurances of m_s.varp or
        // m_s.expr
        visit(specfuncp);
    }

public:
    explicit SpeculativeReplaceVisitor(AstNodeModule* modp, AstMTaskBody* bodyp, Speculateable& s,
                                       bool branch)
        : m_s(s)
        , m_branch(branch)
        , m_mtaskp(bodyp)
        , m_modp(modp) {
        visit(bodyp);
    }

    const VarReplMapping& specOutVars() const { return m_specOutVars; }
    const std::set<AstCFunc*>& replacedFunctions() const { return m_replacedFunctions; }
};

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

//######################################################################
// V3Speculation

V3Speculation::V3Speculation() {
    // Process each module in turn
    for (AstNodeModule* nodep = v3Global.rootp()->modulesp(); nodep;
         nodep = VN_CAST(nodep->nextp(), NodeModule)) {
        speculateModule(nodep);
    }
}

AstNode* genCommitSpecVarStmts(const VarReplMapping& mapping) {
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

void emitSpecResolution(AstMTaskBody* bodyp, ExecMTask* mtaskp, bool branch,
                        const Speculateable& s, VarReplMapping& replacedVars, string hierName) {
    AstNode* condp = nullptr;
    if (s.specBoolVar) {
        condp = new AstVarRef(s.specBoolVar->fileline(), s.specBoolVar, VAccess::READ);
        static_cast<AstVarRef*>(condp)->hiernameToProt(hierName);
    } else if (s.condSpec.condExpr) {
        condp = s.condSpec.condExpr->cloneTree(false);
    } else {
        assert(false);
    }

    if (!branch) { condp = new AstNot(condp->fileline(), condp); }
    bodyp->addStmtsp(new AstSpecResolveBool(s.cons->bodyp()->fileline(), mtaskp, condp,
                                            genCommitSpecVarStmts(replacedVars)));
}

void loopDetect() {
    v3Global.rootp()->execGraphp()->mutableDepGraphp()->dumpDotFilePrefixedAlways("loopcheck");
    // DEBUG: loop detection
    v3Global.rootp()->execGraphp()->mutableDepGraphp()->order();
}

void V3Speculation::doSpeculation(AstNodeModule* modp, Speculateable s) {
    auto* execGraphp = v3Global.rootp()->execGraphp();
    updateDataflowInfo(modp);

    assert(static_cast<bool>(s.specBoolVar) ^ static_cast<bool>(s.condSpec.brExpr)
           && "A speculation must be either boolean speculation or conditional speculation");

    GatherScopeNamesVisitor scopeNames(s.cons);
    string scopeName;
    if (s.specBoolVar) { scopeName = scopeNames.at(s.specBoolVar); }
    m_dfgs[s.cons]->dumpDotFilePrefixedAlways("DFG_" + cvtToStr(s.cons->id()));

    // Before duplicating MTask, perform replacement of conditional expression. This method will
    // replace the conditional expr with a temporary variable. This variable will then be replaced
    // by SpeculativeReplaceVisitor.
    if (s.isCondSpec()) {
        SpecCondExprReplaceVisitor specExprVisitor(s.cons->bodyp(), s);
        s.condSpec.dummyVar = specExprVisitor.dummyVar();
        s.condSpec.condExpr = specExprVisitor.condExpr();
        assert(s.condSpec.condExpr && s.condSpec.dummyVar);
    }

    // Duplicate the MTask, representing the true and false constant value of the speculation
    // variable
    auto* consmtbodyp_t = s.cons->bodyp()->cloneTree(false);
    auto* consmtbodyp_f = s.cons->bodyp()->cloneTree(false);

    // Speculate into them. Manage scope of SpeculativeReplaceVisitors to ensure destruction of
    // user pointers.
    std::set<AstCFunc*> replacedFunctions;
    auto performSpeculation = [&](AstMTaskBody* mtaskbody, bool branch) {
        { AstDotDumper(mtaskbody).dumpDotFilePrefixedAlways("specbody"); }
        SpeculativeReplaceVisitor specVisitor(modp, mtaskbody, s, branch);
        replacedFunctions.insert(specVisitor.replacedFunctions().begin(),
                                 specVisitor.replacedFunctions().end());

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
    consmtp_t->speculative(ExecMTask::Speculative::True, s.cons->id(), consmtp_f);
    consmtp_f->speculative(ExecMTask::Speculative::False, s.cons->id(), consmtp_t);
    consmtbodyp_t->execMTaskp(consmtp_t);
    consmtbodyp_f->execMTaskp(consmtp_f);

    s.prod->addDownstreamSpeculativeTask(consmtp_t);
    s.prod->addDownstreamSpeculativeTask(consmtp_f);

    // Inherit replaced MTasks' downstream tasks
    for (const auto& mtaskp : s.cons->downstreamSpeculativeMTasks()) {
        consmtp_f->addDownstreamSpeculativeTask(mtaskp);
        consmtp_t->addDownstreamSpeculativeTask(mtaskp);
    }

    // Create dependency edges for all incoming variables except the speculated boolean. We cannot
    // use the dependency graph here, due to transitive edge removal, and by the fact that the only
    // remaining dependency edge in the graph is the actual boolean which we are currently
    // speculating on.
    // @todo: Only create edge if actual dependency (ie. avoid loops, we need topological sort)
    for (const auto& invarp : m_io.at(modp).at(s.cons).ins) {
        if (s.specVars.count(invarp) != 0) continue;
        if (invarp->varType() == AstVarType::PORT) continue;

        auto prodMTaskp = m_varProducedBy.find(invarp);
        if (prodMTaskp == m_varProducedBy.end()) continue;
        std::vector<ExecMTask*> inEdges = {prodMTaskp->second};
        if (prodMTaskp->second->partnerSpecMTaskp() != nullptr) {
            // Also depend on the other speculative branch
            inEdges.push_back(prodMTaskp->second->partnerSpecMTaskp());
        }
        for (auto* prodp : inEdges) {
            assert(prodp != consmtp_t);
            new V3GraphEdge(execGraphp->mutableDepGraphp(), prodp, consmtp_t, 1);
            assert(prodp != consmtp_f);
            new V3GraphEdge(execGraphp->mutableDepGraphp(), prodp, consmtp_f, 1);
        }
    }

    // Add speculative resolution statements to end of speculative bodies
    emitSpecResolution(consmtbodyp_t, consmtp_t, true, s, specTOutVars, scopeName);
    emitSpecResolution(consmtbodyp_f, consmtp_f, false, s, specFOutVars, scopeName);

    // The speculated mtasks should inherit the outgoing edges from their origin MTask
    for (V3GraphEdge* edgep = s.cons->outBeginp(); edgep; edgep = edgep->outNextp()) {
        ExecMTask* toExecMTaskp = dynamic_cast<ExecMTask*>(edgep->top());
        assert(consmtp_t != toExecMTaskp);
        new V3GraphEdge(execGraphp->mutableDepGraphp(), consmtp_t, toExecMTaskp, 1);
        assert(consmtp_f != toExecMTaskp);
        new V3GraphEdge(execGraphp->mutableDepGraphp(), consmtp_f, toExecMTaskp, 1);
    }

    // Finally, delete the original MTask and its speculated functions
    s.cons->bodyp()->unlinkFrBack()->deleteTree(), s.cons->bodyp();
    s.cons->unlinkDelete(execGraphp->mutableDepGraphp());
    for (auto& specFunc : replacedFunctions) { specFunc->unlinkFrBack()->deleteTree(); }

    execGraphp->mutableDepGraphp()->dumpDotFilePrefixedAlways("dep_post_spec");
}

// Locate MTasks which depend on this variable:
// - There is an edge between the two MTasks
// - The given variable is the only depending variable produced by the source MTask
//
// Due to transitive edges, the consuming MTask can inherit the incoming dependencies of
// the producing mtask through speculation.
bool V3Speculation::isCriticalVariable(AstVar* varp, ExecMTask* consp) {
    // If variable is available at evaluation time (ie. register/sync. input), it can never be a
    // dritical variable. This is shown by the lack of producing MTasks.
    if (varp->prodMtaskIds().empty()) { return false; }

    ExecMTask* prodp = m_mtaskIdToMTask.at(*varp->prodMtaskIds().begin());

    // Locate MTasks which depend on this variable:
    // - There is an edge between the two MTasks
    // - The given variable is the only depending variable produced by the source MTask
    //
    // Due to transitive edges, the consuming MTask can inherit the incoming dependencies of
    // the producing mtask through speculation.
    for (V3GraphEdge* edgep = prodp->outBeginp(); edgep; edgep = edgep->outNextp()) {
        if (consp != dynamic_cast<ExecMTask*>(edgep->top())) continue;

        std::vector<AstVar*> sharedVars;

        auto& prodOuts = m_dfgs.at(prodp)->outs();
        auto& consIns = m_dfgs.at(consp)->ins();

        assert(std::find(prodOuts.begin(), prodOuts.end(), varp) != prodOuts.end()
               && "varp should be in produced set of producing mtask");

        std::set_intersection(prodOuts.begin(), prodOuts.end(), consIns.begin(), consIns.end(),
                              std::back_inserter(sharedVars));
        return sharedVars.size() == 1;
    }
    return false;
}

void V3Speculation::gatherBoolVarSpecs(AstNodeModule* modp, Speculateables& speculateable) {
    V3Graph* mtasksp = v3Global.rootp()->execGraphp()->mutableDepGraphp();

    // Locate module variables
    std::vector<AstVar*> boolVars1;
    for (AstNode* nodep = modp->stmtsp(); nodep; nodep = nodep->nextp()) {
        if (AstVar* varp = dynamic_cast<AstVar*>(nodep)) {
            if (varp->widthMin() == 1) { boolVars1.push_back(varp); }
        }
    }

    // Dependency filtering
    // Filter to the set of dependencies which result as an in-evaluation dependency
    // rank(producing node) < rank(consuming node)
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

            // For now, also filter tasks which have more than 1 incoming edge in the dependency
            // graph, post transitive edge removal. This is to avoid speculating when we anyways
            // have to wait on another task to finish. @Todo: This should be guided by scheduling!
            int nonTransIn = 0;
            for (auto* edgep = consp->inBeginp(); edgep; edgep = edgep->inNextp()) {
                nonTransIn++;
            }
            if (nonTransIn > 1) { continue; }

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
                s.specBoolVar = varp;
                s.specVars.insert(varp);
                speculateable[s.cons].push_back(s);
            }
        }
    }
}

class ConditionalSpecDetectionVisitor final : public AstNVisitor {
private:
    struct Candidate {
        // Set of variables contained within the possibly speculative conditional
        std::set<AstVar*> specVars;
        AstNode* brExpr = nullptr;
    };

    // A stack of the current candidates that we are traversing
    std::vector<Candidate*> inCandidateStack;

    std::set<AstVar*> reffedNonSpecVars;
    std::vector<Candidate> m_candidates;

    // Set of variables referenced outside the speculated variable
    std::set<AstVar*> refNonSpecVars;

    Candidate const* chosenCandidate = nullptr;

public:
    void visit(AstNode* nodep) { iterateChildren(nodep); }

    void visit(AstVarRef* varrefp) {
        if (!inCandidateStack.empty()) {
            // For all candidates in the stack, this will be a speculative reference
            for (auto& candidate : inCandidateStack) {
                candidate->specVars.insert(varrefp->varp());
            }
        } else {
            reffedNonSpecVars.insert(varrefp->varp());
        }
        visit(static_cast<AstNode*>(varrefp));
    }

    void handleCondExpr(AstNode* nodep, AstNode* ifp, AstNode* thenp, AstNode* elsep) {
        m_candidates.emplace_back();
        auto* newCandidate = &(*m_candidates.rbegin());
        inCandidateStack.push_back(newCandidate);
        newCandidate->brExpr = nodep;
        visit(static_cast<AstNode*>(ifp));
        inCandidateStack.pop_back();
        if (thenp) visit(static_cast<AstNode*>(thenp));
        if (elsep) visit(static_cast<AstNode*>(elsep));
    }

    void visit(AstNodeCond* nodep) {
        handleCondExpr(nodep, nodep->condp(), nodep->expr1p(), nodep->expr2p());
    }

    void visit(AstNodeIf* nodep) {
        handleCondExpr(nodep, nodep->condp(), nodep->ifsp(), nodep->elsesp());
    }

    void visit(AstCCall* ccallp) { visit(ccallp->funcp()); }

public:
    explicit ConditionalSpecDetectionVisitor(AstMTaskBody* bodyp) {
        // Traverse the circuit to detect candidate conditional expressions
        visit(bodyp);
        // Determine the conditional expression where [specVars \ nonSpecRefVars = Ø]
        // The expression should be the earliest in the list of candidates, since they are
        // traversed by depth in the expression tree.
        for (const auto& candidate : m_candidates) {
            std::vector<AstVar*> intersection;
            std::set_intersection(candidate.specVars.begin(), candidate.specVars.end(),
                                  refNonSpecVars.begin(), refNonSpecVars.end(),
                                  std::back_inserter(intersection));
            if (intersection.empty()) {
                chosenCandidate = &candidate;
                break;
            }
        }
    }

    bool isSpeculateable() const { return chosenCandidate != nullptr; }
    std::set<AstVar*> condVars() const {
        return chosenCandidate ? chosenCandidate->specVars : std::set<AstVar*>();
    }
    AstNode* brExpr() const { return chosenCandidate ? chosenCandidate->brExpr : nullptr; }
};

void V3Speculation::gatherConditionalSpecs(AstNodeModule* modp, Speculateables& speculateables) {
    V3Graph* mtasksp = v3Global.rootp()->execGraphp()->mutableDepGraphp();

    // For each mtask, ...
    for (auto* vtxp = mtasksp->verticesBeginp(); vtxp; vtxp = vtxp->verticesNextp()) {
        auto* mtaskp = static_cast<ExecMTask*>(vtxp);
        ConditionalSpecDetectionVisitor v(mtaskp->bodyp());

        // Speculate only if isSpeculateable and it actually has an incoming dependency edge.
        if (v.isSpeculateable() && (mtaskp->inBeginp() != nullptr)) {
            // Speculate if one of the variables inside the speculated conditional expression is a
            // critical variable
            for (auto* varp : v.condVars()) {
                if (isCriticalVariable(varp, mtaskp)) {
                    Speculateable s;
                    s.prod = m_mtaskIdToMTask.at(*varp->prodMtaskIds().begin());
                    s.cons = mtaskp;
                    s.condSpec.brExpr = v.brExpr();
                    s.specVars = v.condVars();
                    speculateables[s.cons].push_back(s);
                    break;
                }
            }
        }
    }
}

void V3Speculation::speculateModule(AstNodeModule* modp) {
    V3Graph* mtasksp = v3Global.rootp()->execGraphp()->mutableDepGraphp();
    updateDataflowInfo(modp);
    std::map<ExecMTask*, std::vector<Speculateable>> speculateable;

    AstDotDumper dump(v3Global.rootp());
    dump.dumpDotFilePrefixedAlways("ast_pre");

    gatherBoolVarSpecs(modp, speculateable);
    gatherConditionalSpecs(modp, speculateable);

    // For now, only speculate a single variable in partitions which could in theory speculate
    // multiple
    for (auto& it : speculateable) {
        if (it.second.size() > 1) { it.second = {it.second.at(0)}; }
    }

    // Go speculate!
    for (const auto& s : speculateable) { doSpeculation(modp, s.second.at(0)); }

    // Run constant propagation to clean up speculate partitions
    V3Const::constifyAll(v3Global.rootp());

    // New mtasks were inserted, reassure order of graph
    mtasksp->order();
}

void V3Speculation::updateDataflowInfo(AstNodeModule* modp) {
    m_dfgs.clear();
    m_io.clear();
    m_mtaskIdToMTask.clear();
    m_varProducedBy.clear();

    // Ensure up-to-date ordering (rank) of MTasks
    v3Global.rootp()->execGraphp()->mutableDepGraphp()->order();

    std::vector<AstVar*> variables;
    for (AstNode* nodep = modp->stmtsp(); nodep; nodep = nodep->nextp()) {
        if (AstVar* varp = dynamic_cast<AstVar*>(nodep)) {
            // Reset producer/consumer IDs (updated later)
            varp->clearMTaskIds();
            variables.push_back(varp);
        }
    }

    // Gather DFGs for each MTask. This also updates the producer/consumer info of each
    // variable
    V3Graph* mtasksp = v3Global.rootp()->execGraphp()->mutableDepGraphp();
    for (V3GraphVertex* vxp = mtasksp->verticesBeginp(); vxp; vxp = vxp->verticesNextp()) {
        ExecMTask* mtaskp = dynamic_cast<ExecMTask*>(vxp);
        m_nextMTaskID = m_nextMTaskID <= mtaskp->id() ? mtaskp->id() + 1 : m_nextMTaskID;
        m_dfgs[mtaskp] = new DFG(mtaskp->bodyp());
        m_mtaskIdToMTask[mtaskp->id()] = mtaskp;
    }

    for (const auto& varp : variables) {
        if (varp->prodMtaskIds().size() == 0) {
            // Probably an input...
            continue;
        }

        // Collect inputs; An input is a dependent variable which is produced by an MTask that has
        // a lower rank than the consuming mtask.
        // If multiple producers, the variable is an input to an MTask if there exists a producing
        // mtask which has lower rank than the consuming mtask.
        for (const int consMTaskid : varp->consMtaskIds()) {
            auto consMTask = m_mtaskIdToMTask.find(consMTaskid);
            assert(consMTask != m_mtaskIdToMTask.end());

            for (const auto& prodMTaskID : varp->prodMtaskIds()) {
                auto prodMTask = m_mtaskIdToMTask.find(prodMTaskID);
                if (prodMTask->second->rank() < consMTask->second->rank()) {
                    m_io[modp][consMTask->second].ins.insert(varp);
                    break;
                }
            }
        }

        // Collect outputs
        for (const int prodMTaskid : varp->prodMtaskIds()) {
            auto it = m_mtaskIdToMTask.find(prodMTaskid);
            assert(it != m_mtaskIdToMTask.end());
            ExecMTask* prodp = it->second;
            m_io[modp][prodp].outs.insert(varp);

            auto* mtask = m_mtaskIdToMTask.at(prodMTaskid);
            assert(mtask->id() == prodMTaskid
                   && "Mismatched mtask ID and actual produced mtask ID");
            m_varProducedBy[varp] = mtask;
        }
    }

    /*
    cout << "IO is:" << endl;
    for (const auto& it : m_io.at(modp)) {
        cout << it.first->name();
        cout << "\t ins: ";
        for (const auto& it2 : it.second.ins) { cout << it2->name(); }
        cout << endl << "\t outs: ";
        for (const auto& it2 : it.second.outs) { cout << it2->name(); }
        cout << endl;
    }
    */
}
