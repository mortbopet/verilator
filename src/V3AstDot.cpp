// -*- mode: C++; c-file-style: "cc-mode" -*-
//*************************************************************************
// DESCRIPTION: Verilator: Dot file dumping of ast nodes
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

#include "V3AstDot.h"
#include "V3File.h"
#include "V3Graph.h"

#include <memory>

//======================================================================

AstDotDumper::AstDotDumper(AstNode* topp)
    : m_topp(topp) {}

void AstDotDumper::dumpDotFile(const std::string& filename) {
    auto logusp = std::unique_ptr<std::ofstream>(V3File::new_ofstream(filename));
    m_logp = logusp.get();
    if (m_logp->fail()) v3fatal("Can't write " << filename);

    // Header
    *m_logp << "digraph v3graph {\n";
    *m_logp << "  graph[labelloc=t labeljust=l label=\"" << filename << "\"]\n";

    iterate(m_topp);

    // Trailer
    *m_logp << "}\n";
    m_logp->close();
}

void AstDotDumper::dumpDotFilePrefixed(const string& nameComment) {
    if (v3Global.opt.dumpTree()) {
        dumpDotFile(v3Global.debugFilename("mt_liveness_" + nameComment) + ".dot");
    }
}

//! Variant of dumpDotFilePrefixed without --dump option check
void AstDotDumper::dumpDotFilePrefixedAlways(const string& nameComment) {
    dumpDotFile(v3Global.debugFilename("mt_liveness_" + nameComment) + ".dot");
}

string AstDotDumper::astNodeDotLabel(AstNode* nodep) {
    auto str = std::ostringstream();
    nodep->dump(str);
    auto quotedStr = VString::quoteStringLiteralForShell(str.str());
    for (size_t i = 0; i < quotedStr.size(); i++) {
        if (quotedStr.at(i) == ' ') quotedStr[i] = '\n';
    }

    return quotedStr;
}

void AstDotDumper::dumpNode(AstNode* nodep) {
    if (m_dumpedIDs.count(nodep) != 0) { return; }

    *m_logp << "  " << astNodeDotName(nodep) << " [label=" << astNodeDotLabel(nodep) << "]\n";
    m_dumpedIDs.insert(nodep);
}

void AstDotDumper::visit(AstNode* nodep) {
    dumpNode(nodep);
    AstDotEdgeVisitor edgevisit(m_logp, nodep);
    edgevisit.iterateChildren(nodep);
    iterateChildren(nodep);
}

void AstDotDumper::visit(AstCCall* nodep) {
    visit(static_cast<AstNode*>(nodep));

    if (followCCalls) {
        AstDotEdgeVisitor edgevisit(m_logp, nodep);
        edgevisit.visit(nodep->funcp());  // Create edge between ccall and cfunc
        visit(nodep->funcp());  // follow c-func
    }
}

AstDotDumper::AstDotEdgeVisitor::AstDotEdgeVisitor(std::ofstream* logp, AstNode* fromp)
    : m_logp(logp)
    , m_fromp(fromp) {}

void AstDotDumper::AstDotEdgeVisitor::visit(AstNode* nodep) {
    *m_logp << "  " << astNodeDotName(m_fromp) << "->" << astNodeDotName(nodep);
    if (nodep->hasDType()) {
        *m_logp << " [label=\"\\[" << nodep->widthMin() - 1 << ":" << 0 << "\\]\"]";
    }
    *m_logp << "\n";
}

AstToDfg::AstToDfg(AstNode* topp)
    : m_topp(topp) {
    iterate(topp);
}

void AstToDfg::visit(AstNode* nodep) {
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

        if (srcvtp != nodeDfgp) { new DFGEdge(this, srcvtp, nodeDfgp, 1); }
    }
}

void AstToDfg::visit(AstVar* nodep) {
    UASSERT(false, "???");

    visit(static_cast<AstNode*>(nodep));
}

void AstToDfg::visit(AstVarRef* nodep) {
    UASSERT(true, "");
    visit(static_cast<AstNode*>(nodep));
}

DFGVertex* AstToDfg::updateVarDFGSrc(AstNode* sourcep, AstVar* varp) {
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

void AstToDfg::visit(AstNodeAssign* nodep) {
    // Ensure depth-first traversal before updating the variable assignment node. As a reminder,
    // generating the DFG from the AST is in general inverting edge directions of the AST dag.
    visit(static_cast<AstNode*>(nodep));

    AstNodeVarRef* varrefp = dynamic_cast<AstNodeVarRef*>(nodep->lhsp());
    UASSERT(varrefp, "Assignment to non-variable?");
    AstVar* varp = varrefp->varp();

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

void AstToDfg::visit(AstCCall* nodep) {
    visit(static_cast<AstNode*>(nodep));
    visit(nodep->funcp());  // follow c-func
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
