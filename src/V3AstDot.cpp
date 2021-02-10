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
    { AstDotEdgeVisitor edgevisit(m_logp, nodep); }
    iterateChildren(nodep);
}

AstDotDumper::AstDotEdgeVisitor::AstDotEdgeVisitor(std::ofstream* logp, AstNode* fromp)
    : m_logp(logp)
    , m_fromp(fromp) {
    iterateChildren(fromp);
}

void AstDotDumper::AstDotEdgeVisitor::visit(AstNode* nodep) {
    *m_logp << "  " << astNodeDotName(m_fromp) << "->" << astNodeDotName(nodep);
    if (nodep->hasDType()) {
        *m_logp << " [label=\"\\[" << nodep->widthMin() - 1 << ":" << 0 << "\\]\"]";
    }
    *m_logp << "\n";
}
