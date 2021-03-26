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
#include "V3Dataflow.h"

struct VarReplInfo {
    // Maintain point of assignment for cloning assignment attributes (fileline etc.)
    AstNodeAssign* assignp = nullptr;

    // Maintain point of reference for cloning attributes (scope etc.)
    AstVarRef* varrefp = nullptr;

    // Replacement variable
    AstVar* replvarp = nullptr;
};

using VarReplMapping = std::unordered_map<AstVar*, VarReplInfo>;

struct CondSpeculateable {
    // Branch expression identified to be speculateable
    AstNode* brExpr = nullptr;
    // Dummy var generated when exchanging old cond expression
    AstVar* dummyVar = nullptr;
    // Old cond expression, that is reemitted during speculative commit
    AstNode* condExpr = nullptr;
};

/**
 * @brief The Speculateable struct
 * Maintains the set of producer/consumer MTasks which share a boolean variable, wherein the
 * relationship has been determined to be speculateable.
 */
struct Speculateable {
    // Set of variables that is referenced within this speculation (which, in turn, are
    // speculative and must be tracked)
    std::set<AstVar*> specVars;

    // A boolean variable that is speculative
    AstVar* specBoolVar = nullptr;
    CondSpeculateable condSpec;
    ExecMTask* prod = nullptr;
    ExecMTask* cons = nullptr;
    bool isCondSpec() const { return condSpec.brExpr != nullptr; }
    bool isBoolSpec() const { return condSpec.brExpr == nullptr; }
};
using Speculateables = std::map<ExecMTask*, std::vector<Speculateable>>;
class V3Speculation {
public:
    V3Speculation();

private:
    struct VarIO {
        std::set<AstVar*> ins;
        std::set<AstVar*> outs;
    };

    AstNode* genCommitSpecVarStmts(const VarReplMapping&);
    void updateDataflowInfo(AstNodeModule* modp);
    void speculateModule(AstNodeModule* nodep);
    void doSpeculation(AstNodeModule* modp, Speculateable s);
    void gatherBoolVarSpecs(AstNodeModule* modp, Speculateables& speculateables);
    void gatherConditionalSpecs(AstNodeModule* modp, Speculateables& speculateables);
    bool isCriticalVariable(AstVar* varp, ExecMTask* consp);

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
    std::unordered_map<AstNodeModule*, std::unordered_map<const ExecMTask*, VarIO>> m_io;
    std::map<AstVar*, ExecMTask*> m_varProducedBy;
    std::map<int, ExecMTask*> m_mtaskIdToMTask;
    unsigned m_nextMTaskID = 0;
};

#endif  // Guard
