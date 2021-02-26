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

#ifndef _V3SPECULATION_H_
#define _V3SPECULATION_H_

#include "V3Ast.h"
#include "V3Graph.h"

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

/**
 * @todo:
 * - maintain modpointer during visitor; emit speculative variables in current module
 * -
 */

class V3Speculation {
public:
    V3Speculation();

private:
    /**
     * @brief The Speculateable struct
     * Maintains the set of producer/consumer MTasks which share a boolean variable, wherein the
     * relationship has been determined to be speculateable.
     */
    struct Speculateable {
        AstVar* specVar = nullptr;
        ExecMTask* prod = nullptr;
        std::vector<ExecMTask*> cons;
    };

    void speculateModule(AstNodeModule* nodep);
    void doSpeculation(AstNodeModule* modp, const Speculateable& s);

    /**
     * @brief removeDependency
     * Removes dependency on @param varp from @param mtaskp. In doing so, we re-create dependency
     * edges from _all other_ variables that mtaskp is dependent on. Later, through removal of
     * transitive edges, only a single of these edges will remain.
     *
     */
    void removeDependency(ExecMTask* mtaskp, AstVar* varp);

    /**
     * @brief m_dataflow
     * Maintains the in- and output variables of MTasks and CFuncs
     */
    std::unordered_map<ExecMTask*, DFG*> m_dfgs;
    std::map<int, ExecMTask*> m_mtaskIdToMTask;
    unsigned m_nextMTaskID = 0;
};

#endif  // Guard
