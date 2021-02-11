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

class DFGVertex VL_NOT_FINAL : public V3GraphVertex {
public:
    DFGVertex(V3Graph* graph, AstNode* nodep)
        : V3GraphVertex(graph)
        , m_nodep(nodep) {}
    virtual ~DFGVertex() override = default;

    virtual string dotName() const override { return v3Global.ptrToId(m_nodep); };
    virtual string name() const override;
    virtual string dotShape() const override;

private:
    AstNode* m_nodep;
};

class DFGEdge VL_NOT_FINAL : public V3GraphEdge {
public:
    DFGEdge(V3Graph* graph, DFGVertex* fromp, DFGVertex* top, int weight)
        : V3GraphEdge(graph, fromp, top, weight) {}
    virtual ~DFGEdge() override = default;

    virtual string dotColor() const override { return color; }
    virtual string dotStyle() const override {
        return style == "" ? V3GraphEdge::dotStyle() : style;
    }

    string style = "";
    string color = "black";

private:
    AstNode* m_nodep;
};

/**
 * @brief The AstToDfg class
 *
 * Strong assumption on AST graph being traversed in order of dependency (ie. CFunc statements
 * traversed in order)
 */

class AstToDfg final : public AstNVisitor, public V3Graph {
public:
    // CONSTRUCTORS
    explicit AstToDfg(AstNode* nodep);

    // Configuration
    bool followCCalls = true;

private:
    virtual void visit(AstNode* nodep) override;
    virtual void visit(AstVar* nodep) override;
    virtual void visit(AstVarRef* nodep) override;
    virtual void visit(AstNodeAssign* nodep) override;
    virtual void visit(AstCCall* nodep) override;
    /**
     * @brief updateVarDFGSrc
     * @param sourcep
     *      Pointer to AstNode which updated the variable reference. If nullptr, varp is
     *      interpreted as being an input variable to the mtask.
     * @param varp
     * @return
     */
    DFGVertex* updateVarDFGSrc(AstNode* sourcep, AstVar* varp);

    AstNode* m_topp;

    // Tracking of current variable to most-recent source node in DFG
    std::map<AstVar*, DFGVertex*> m_varToDFGVp;
    std::map<AstNode*, DFGVertex*> m_nodeToDFGVp;
};

#endif  // Guard
