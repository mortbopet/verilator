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

#include <algorithm>

//######################################################################
// Utilities

class SpecVarMarkVisitor final : public AstNVisitor {
    // NODE STATE
    // AstVar::user4p           -> bool, Var marked, 0=not set yet
private:
    // VISITORS
    virtual void visit(AstVarRef* nodep) override {
        if (nodep->varp()) nodep->varp()->user4(1);
    }
    virtual void visit(AstNode* nodep) override { iterateChildren(nodep); }

public:
    // CONSTRUCTORS
    explicit SpecVarMarkVisitor(AstNode* nodep) {
        AstNode::user4ClearTree();  // Check marked InUse before we're called
        iterate(nodep);
    }
    virtual ~SpecVarMarkVisitor() override = default;
};

class SpecVarFindVisitor final : public AstNVisitor {
    // NODE STATE
    // AstVar::user4p           -> bool, input from SpecVarMarkVisitor
    // MEMBERS
    bool m_found = false;

private:
    // VISITORS
    virtual void visit(AstVarRef* nodep) override {
        if (nodep->varp() && nodep->varp()->user4()) m_found = true;
    }
    virtual void visit(AstNode* nodep) override { iterateChildren(nodep); }

public:
    // CONSTRUCTORS
    explicit SpecVarFindVisitor(AstNode* nodep) { iterateAndNextNull(nodep); }
    virtual ~SpecVarFindVisitor() override = default;
    // METHODS
    bool found() const { return m_found; }
};

//######################################################################
// DFG

DFG::DFG(AstMTaskBody* bodyp)
    : m_bodyp(bodyp) {
    iterate(bodyp);
}

void DFG::visit(AstNode* nodep) {
    if (!nodep) return;

    iterateChildren(nodep);

    DFGVertex* nodeDfgp = nullptr;

    for (auto* childp : {nodep->op1p(), nodep->op2p(), nodep->op3p(), nodep->op4p()}) {
        if (childp == nullptr) continue;
        DFGVertex* srcvtp = nullptr;
        auto* varrefp = dynamic_cast<AstVarRef*>(childp);
        if (varrefp &&
            /* handled in visit(AstNodeAssign) */
            dynamic_cast<AstNodeAssign*>(nodep) == nullptr) {
            // Locate latest write to varref
            if (m_varToDFGVp.count(varrefp->varp()) == 0) {
                // First time seeing variable
                updateVarDFGSrc(nullptr, varrefp->varp());
            }
            srcvtp = m_varToDFGVp.at(varrefp->varp());
        } else {
            if (m_nodeToDFGVp.count(nodep) == 0) {
                m_nodeToDFGVp[nodep] = new DFGVertex(this, nodep);
            }
            nodeDfgp = m_nodeToDFGVp[nodep];
            if (m_nodeToDFGVp.count(childp) == 0) {
                m_nodeToDFGVp[childp] = new DFGVertex(this, childp);
            }
            srcvtp = m_nodeToDFGVp[childp];
        }

        if (varrefp) {
            m_io.ins.insert(varrefp->varp());
            varrefp->varp()->addConsumingMTaskId(m_bodyp->execMTaskp()->id());
        }
        if (srcvtp != nodeDfgp && nodeDfgp != nullptr) { new DFGEdge(this, srcvtp, nodeDfgp, 1); }
    }
}

void DFG::visit(AstVar* nodep) {
    // UASSERT(false, "???");
    visit(static_cast<AstNode*>(nodep));
}

void DFG::visit(AstVarRef* nodep) {
    UASSERT(true, "");
    visit(static_cast<AstNode*>(nodep));
}

DFGVertex* DFG::updateVarDFGSrc(AstNode* sourcep, AstVar* varp) {
    DFGVertex* varSrcp = nullptr;

    auto dfgvp = m_varToDFGVp.find(varp);
    if (dfgvp == m_varToDFGVp.end()) {
        // First reference to variable
        varSrcp = new DFGVertex(this, varp);
        m_nodeToDFGVp[varp] = varSrcp;
    } else {
        UASSERT(sourcep != nullptr, "Should always be present");
        varSrcp = m_nodeToDFGVp.at(sourcep);
        // varSrcp = new DFGVertex(this, sourcep);
        // m_nodeToDFGVp[sourcep] = varSrcp;
    }

    m_varToDFGVp[varp] = varSrcp;
    return varSrcp;
}

void DFG::visit(AstNodeAssign* nodep) {
    // Ensure depth-first traversal before updating the variable assignment node. As a reminder,
    // generating the DFG from the AST is in general inverting edge directions of the AST dag.
    visit(static_cast<AstNode*>(nodep));

    AstNodeVarRef* varrefp = dynamic_cast<AstNodeVarRef*>(nodep->lhsp());
    if (varrefp == nullptr) { return; }
    UASSERT(varrefp, "Assignment to non-variable?");

    AstVar* varp = varrefp->varp();

    m_io.outs.insert(varp);
    varrefp->varp()->addProducingMTaskId(m_bodyp->execMTaskp()->id());

    auto it = m_varToDFGVp.find(varp);
    DFGVertex* oldsrcp = nullptr;

    if (it == m_varToDFGVp.end()) {
        // First time seeing variable
        oldsrcp = updateVarDFGSrc(nullptr, varrefp->varp());
    } else {
        oldsrcp = it->second;
    }
    // This is now the most-recent write to the variable
    auto thisSrcp = updateVarDFGSrc(nodep, varp);

    // Create edge from old src to this assignment
    auto varRefEdge = new DFGEdge(this, oldsrcp, thisSrcp, 1);
    varRefEdge->style = "dashed";
    varRefEdge->color = "red";
}

void DFG::visit(AstCCall* nodep) {
    visit(static_cast<AstNode*>(nodep));
    visit(nodep->funcp());  // follow c-func
}

void DFG::visit(AstNodeIf* nodep) {
    visit(static_cast<AstNode*>(nodep));

    // In order, handle condition, then-branch and else-branch
    // visit(nodep->condp());
    // visit(nodep->ifsp());
    // visit(nodep->elsesp());
}

string DFGVertex::dotShape() const {
    if (dynamic_cast<AstConst*>(m_nodep)) {
        return "square";
    } else if (auto* varp = dynamic_cast<AstVar*>(m_nodep)) {
        return "diamond";
    }

    return "";
}

string DFGVertex::name() const {
    if (auto* assp = dynamic_cast<AstNodeAssign*>(m_nodep)) {
        if (auto* varp = dynamic_cast<AstVarRef*>(assp->lhsp())) {
            return varp->varp()->origName() + " <=";
        } else {
            return "<=";
        }
    } else if (auto* varp = dynamic_cast<AstVar*>(m_nodep)) {
        return varp->origName();
    } else if (auto* constp = dynamic_cast<AstConst*>(m_nodep)) {
        return constp->name();
    } else if (m_nodep == nullptr) {
        return "sink";
    } else {
        return m_nodep->typeName();
    }
}

//######################################################################
// SpeculationReplaceVisitor

inline string specNameSuffix(bool branch) { return "__SPEC__" + string(branch ? "t" : "f"); }

class SpeculativeReplaceVisitor final : public AstNVisitor {
private:
    // Mapping between non-speculative to replaced speculative variable, which is being written to
    // in the partition.
    std::unordered_map<AstVar*, AstVar*> m_specOutVars;
    AstVar* m_varp;  // Boolean variable to speculate on
    bool m_branch;
    AstMTaskBody* m_mtaskp;
    AstNodeModule* m_modp;
    int m_it = 0;

    string nameSuffix() const { return specNameSuffix(m_branch); }

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

    void visit(AstNode* nodep) { iterateChildren(nodep); }

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
        AstDotDumper dumper(m_mtaskp);
        dumper.dumpDotFilePrefixedAlways("pre_ccall" + std::to_string(m_it));

        // Create a duplicate of the called C function, and perform boolean const propagation
        // into this function.
        const AstCFunc* origfuncp = ccallp->funcp();
        AstCFunc* specfuncp = ccallp->funcp()->cloneTree(true);  // i think true here?
        specfuncp->name(specfuncp->name() + nameSuffix());
        origfuncp->scopep()->addActivep(specfuncp);

        // Replace the function call to the new, speculative function
        auto newFuncCall = new AstCCall(ccallp->fileline(), specfuncp, ccallp->argsp());
        ccallp->replaceWith(newFuncCall);

        AstDotDumper dumper2(m_mtaskp);
        dumper2.dumpDotFilePrefixedAlways("post_ccall" + std::to_string(m_it));

        m_it++;

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
        // Duplicate the MTask, representing the true and false constant value of the speculation
        // variable
        auto* consmtbodyp_t = consmtaskp->bodyp()->cloneTree(false);
        auto* consmtbodyp_f = consmtaskp->bodyp()->cloneTree(false);

        // Speculate into them
        SpeculativeReplaceVisitor(modp, consmtbodyp_t, s.specVar, true);
        SpeculativeReplaceVisitor(modp, consmtbodyp_f, s.specVar, false);

        // Insert them into the dependency graph
        execGraphp->addMTaskBody(consmtbodyp_t);
        execGraphp->addMTaskBody(consmtbodyp_f);

        // Create ExecMTasks
        auto* consmtp_t = new ExecMTask(execGraphp->mutableDepGraphp(), consmtbodyp_t,
                                        m_nextMTaskID++, ExecMTask::Speculative::True);
        auto* consmtp_f = new ExecMTask(execGraphp->mutableDepGraphp(), consmtbodyp_f,
                                        m_nextMTaskID++, ExecMTask::Speculative::False);

        // Create dependency edges for all incoming variables except the speculated boolean
    }

    // Remove redundant dependency edges
    execGraphp->mutableDepGraphp()->removeTransitiveEdges();
}

void V3Speculation::speculateModule(AstNodeModule* modp) {
    // Locate module variables
    std::vector<AstVar*> variables;
    for (AstNode* nodep = modp->stmtsp(); nodep; nodep = nodep->nextp()) {
        if (AstVar* varp = dynamic_cast<AstVar*>(nodep)) {
            // Reset producer/consumer IDs (updated later)
            varp->clearMTaskIds();
            variables.push_back(varp);
        }
    }

    // Gather DFGs for each MTask. This also updates the producer/consumer info of each variable
    const V3Graph* mtasksp = v3Global.rootp()->execGraphp()->depGraphp();
    for (V3GraphVertex* vxp = mtasksp->verticesBeginp(); vxp; vxp = vxp->verticesNextp()) {
        ExecMTask* mtaskp = dynamic_cast<ExecMTask*>(vxp);
        m_nextMTaskID = m_nextMTaskID <= mtaskp->id() ? mtaskp->id() + 1 : m_nextMTaskID;
        m_dfgs[mtaskp] = new DFG(mtaskp->bodyp());
        m_mtaskIdToMTask[mtaskp->id()] = mtaskp;
    }

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
    v3Global.rootp()->execGraphp()->mutableDepGraphp()->order();
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

    return;
}
