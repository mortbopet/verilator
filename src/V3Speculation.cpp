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

DFG::DFG(AstNode* topp)
    : m_topp(topp) {
    iterate(topp);
}

void DFG::visit(AstNode* nodep) {
    if (!nodep) return;

    if (m_nodeToDFGVp.count(nodep) == 0) { m_nodeToDFGVp[nodep] = new DFGVertex(this, nodep); }
    auto* nodeDfgp = m_nodeToDFGVp.at(nodep);
    iterateChildren(nodep);

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
            auto it = m_nodeToDFGVp.find(childp);
            UASSERT(it != m_nodeToDFGVp.end(), "Should have been created");
            srcvtp = it->second;
        }

        if (varrefp) { m_io.ins.insert(varrefp->varp()); }

        if (srcvtp != nodeDfgp) { new DFGEdge(this, srcvtp, nodeDfgp, 1); }
    }
}

void DFG::visit(AstVar* nodep) {
    UASSERT(false, "???");

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
    UASSERT(varrefp, "Assignment to non-variable?");
    AstVar* varp = varrefp->varp();

    m_io.outs.insert(varp);

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
    visit(nodep->condp());
    visit(nodep->ifsp());
    visit(nodep->elsesp());
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
// V3Speculation

V3Speculation::V3Speculation() {}

void V3Speculation::go() {
    // Gather DFGs for each MTask
    const V3Graph* mtasksp = v3Global.rootp()->execGraphp()->depGraphp();
    for (V3GraphVertex* vxp = mtasksp->verticesBeginp(); vxp; vxp = vxp->verticesNextp()) {
        ExecMTask* mtaskp = dynamic_cast<ExecMTask*>(vxp);
        m_dfgs[mtaskp] = new DFG(mtaskp->bodyp());
    }

    if (debug()) {
        for (const auto& it : m_dfgs) { it.second->dumpDotFilePrefixedAlways(it.first->name()); }
    }
}
