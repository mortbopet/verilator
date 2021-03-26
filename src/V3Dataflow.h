// -*- mode: C++; c-file-style: "cc-mode" -*-
//*************************************************************************
// DESCRIPTION: Verilator: Speculative partition detection and emission
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

#ifndef _V3DATAFLOW_H_
#define _V3DATAFLOW_H_

#include "V3Ast.h"
#include "V3Graph.h"

struct VarMTaskScope {
    std::set<ExecMTask*> producers;
    std::set<ExecMTask*> consumers;
};

struct MTaskVarIO {
    std::set<AstVar*> ins;
    std::set<AstVar*> outs;
};

class DFGVertex VL_NOT_FINAL : public V3GraphVertex {
public:
    DFGVertex(V3Graph* graph, AstNode* nodep)
        : V3GraphVertex(graph)
        , m_nodep(nodep) {}
    virtual ~DFGVertex() override = default;

    virtual string dotName() const override { return v3Global.ptrToId(m_nodep, false); };
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

class DFG final : public AstNVisitor, public V3Graph {
public:
    // CONSTRUCTORS
    explicit DFG(AstMTaskBody* bodyp);
    DFG() { UASSERT(false, "Illegal null-construction of DFG"); }

    // Configuration
    bool followCCalls = true;

    struct VarIO {
        std::set<AstVar*> ins;
        std::set<AstVar*> outs;
    };

    const std::set<AstVar*>& ins() const { return m_io.ins; }
    const std::set<AstVar*>& outs() const { return m_io.outs; }

private:
    virtual void visit(AstNode* nodep) override;
    virtual void visit(AstVar* nodep) override;
    virtual void visit(AstVarRef* nodep) override;
    virtual void visit(AstNodeAssign* nodep) override;
    virtual void visit(AstCCall* nodep) override;
    virtual void visit(AstNodeIf* nodep) override;

    /**
     * @brief updateVarDFGSrc
     * @param sourcep
     *      Pointer to AstNode which updated the variable reference. If nullptr, varp is
     *      interpreted as being an input variable to the mtask.
     * @param varp
     * @return
     */
    DFGVertex* updateVarDFGSrc(AstNode* sourcep, AstVar* varp);

    AstMTaskBody* m_bodyp;
    VarIO m_io;

    // Tracking of current variable to most-recent source node in DFG
    std::map<AstVar*, DFGVertex*> m_varToDFGVp;
    std::map<AstNode*, DFGVertex*> m_nodeToDFGVp;
};

class V3Dataflow {
public:
    V3Dataflow();

    /**
     * @brief scopeForVar
     * @returns producing and consuming MTasks for @param varp, or nullptr if no scope was
     * recorded.
     */
    const VarMTaskScope* scopeForVar(const AstVar* varp);

    /**
     * @brief ioForMTask
     * @returns in- and output variables for @param mtaskp, or nullptr if no io was recorded.
     */
    const MTaskVarIO* ioForMTask(const ExecMTask* mtaskp);

private:
};

#endif  // Guard