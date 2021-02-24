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

#ifndef _V3ASTDOT_H_
#define _V3ASTDOTH_ 1

#include "V3Ast.h"
#include "V3Global.h"
#include "V3Graph.h"

#include <memory>

//######################################################################

class AstDotDumper final : public AstNVisitor {
public:
    // CONSTRUCTORS
    explicit AstDotDumper(AstNode* nodep);

    // Non-const due to use of AstNVisitor
    void dumpDotFile(const std::string& filename);
    void dumpDotFilePrefixed(const string& nameComment);
    void dumpDotFilePrefixedAlways(const string& nameComment);

    // Configuration
    bool followCCalls = true;

private:
    static string astNodeDotName(AstNode* nodep) { return "n_" + v3Global.ptrToId(nodep, false); }
    string astNodeDotLabel(AstNode* nodep);

    class AstDotEdgeVisitor final : public AstNVisitor {
    public:
        explicit AstDotEdgeVisitor(std::ofstream* logp, AstNode* fromp);
        virtual void visit(AstNode* nodep) override;

    private:
        std::ofstream* m_logp = nullptr;
        AstNode* m_fromp = nullptr;
    };

    virtual void visit(AstNode* nodep) override;
    virtual void visit(AstCCall* nodep) override;

    void dumpNode(AstNode* nodep);

    std::set<AstNode*> m_dumpedIDs;
    std::ofstream* m_logp;
    AstNode* m_topp;
};

#endif  // Guard
