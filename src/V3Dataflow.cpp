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

#include "V3Dataflow.h"
#include "V3Ast.h"
#include "V3AstDot.h"
#include "V3AstNodes.h"
#include "V3PartitionGraph.h"

#include <algorithm>

//######################################################################
// DFG

DFG::DFG(AstMTaskBody* bodyp)
    : m_bodyp(bodyp) {
    iterate(bodyp);

    if (debug() || true) {
        std::cout << "DFG:" << bodyp->execMTaskp()->name() + " I/O is" << endl;
        std::cout << "\tins: ";
        for (const auto& in : m_io.ins) { std::cout << in->origName() << ", "; }
        std::cout << std::endl << "\touts: ";
        for (const auto& in : m_io.ins) { std::cout << in->origName() << ", "; }
        std::cout << std::endl;
    }
}

void DFG::visit(AstNode* nodep) {
    if (!nodep) return;
    auto nodeDfgp = m_nodeToDFGVp.find(nodep);
    if (nodeDfgp == m_nodeToDFGVp.end()) {
        nodeDfgp = m_nodeToDFGVp.insert(m_nodeToDFGVp.end(), {nodep, new DFGVertex(this, nodep)});
    }

    iterateChildren(nodep);
    for (auto* childp : {nodep->op1p(), nodep->op2p(), nodep->op3p(), nodep->op4p()}) {
        if (childp == nullptr) continue;
        DFGVertex* srcvtp = m_nodeToDFGVp[childp];
        UASSERT(srcvtp != nodeDfgp->second, "?");
        new DFGEdge(this, srcvtp, nodeDfgp->second, 1);
    }
}

void DFG::visit(AstVar* nodep) {
    // UASSERT(false, "???");
    visit(static_cast<AstNode*>(nodep));
}

void DFG::visit(AstVarRef* varrefp) {
    visit(static_cast<AstNode*>(varrefp));  // depth-first traversal

    // Locate latest write to varref
    if (m_varToDFGVp.count(varrefp->varp()) == 0) {
        // First time seeing variable
        updateVarDFGSrc(nullptr, varrefp->varp());
    }

    // Register the variable as an input to this MTask if the MTask has not previously
    // written to the variable
    auto* varp = varrefp->varp();
    if (m_io.outs.count(varp) == 0) {
        m_io.ins.insert(varrefp->varp());
        varrefp->varp()->addConsumingMTaskId(m_bodyp->execMTaskp()->id());
    }
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

V3Dataflow::V3Dataflow() { v3Global.rootp()->execGraphp()->mutableDepGraphp(); }
