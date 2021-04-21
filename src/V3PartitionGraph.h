// -*- mode: C++; c-file-style: "cc-mode" -*-
//*************************************************************************
// DESCRIPTION: Verilator: Threading's graph structures
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

#ifndef _V3PARTITIONGRAPH_H_
#define _V3PARTITIONGRAPH_H_

#include "config_build.h"
#include "verilatedos.h"

#include "V3Graph.h"
#include "V3OrderGraph.h"

#include <list>
#include <optional>

//*************************************************************************
// MTasks and graph structures

class AbstractMTask VL_NOT_FINAL : public V3GraphVertex {
public:
    AbstractMTask(V3Graph* graphp)
        : V3GraphVertex{graphp} {}
    virtual ~AbstractMTask() override = default;
    virtual uint32_t id() const = 0;
    virtual uint32_t cost() const = 0;
};

class AbstractLogicMTask VL_NOT_FINAL : public AbstractMTask {
public:
    // TYPES
    typedef std::list<MTaskMoveVertex*> VxList;
    // CONSTRUCTORS
    AbstractLogicMTask(V3Graph* graphp)
        : AbstractMTask{graphp} {}
    virtual ~AbstractLogicMTask() override = default;
    // METHODS
    // Set of logic vertices in this mtask. Order is not significant.
    virtual const VxList* vertexListp() const = 0;
    virtual uint32_t id() const override = 0;  // Unique id of this mtask.
    virtual uint32_t cost() const override = 0;
};

class ExecMTask final : public AbstractMTask {
public:
    enum class Speculative { None, True, False };

private:
    AstMTaskBody* m_bodyp;  // Task body
    uint32_t m_id;  // Unique id of this mtask.
    uint32_t m_priority = 0;  // Predicted critical path from the start of
    // this mtask to the ends of the graph that are reachable from this
    // mtask. In abstract time units.
    uint32_t m_downstreamCost = 0;  // Predicted sum of costs of all descendants
    // of this mtask.
    uint32_t m_cost = 0;  // Predicted runtime of this mtask, in the same
    // abstract time units as priority().
    uint32_t m_startTime = 0;  // predicted start time of this mtask, in the same
    // abstract time units as priority()
    uint32_t m_endTime = 0;  // predicted end time of this mtask, in the same
    // abstract time units as priority()
    uint32_t m_thread = 0xffffffff;  // Thread for static (pack_mtasks) scheduling,
    // or 0xffffffff if not yet assigned.
    const ExecMTask* m_packNextp = nullptr;  // Next for static (pack_mtasks) scheduling
    bool m_threadRoot = false;  // Is root thread
    std::optional<bool> m_spec;  // Speculated branch, if any
    unsigned m_speculatesMtaskId = -1;  // ID of original Mtask that this mtask speculates (we
                                        // don't keep a pointer since the MTask will be removed)

    // If this metask is speculative, the following variable contains the IDs of the MTasks that
    // produces the variable which this MTask speculates on. This can be multiple MTasks, given
    // that it may speculate on a variable produced by other speculative MTasks, which in turn mean
    // that it has a dependency on both of those MTasks.
    std::set<uint32_t> m_upstreamSpeculativeDepMTasks;

    // If this MTask produces a variable which is speculated upon by other mtasks, these
    // speculative mtasks are maintained in the following set. We maintain as IDs rather than
    // pointers, because MTasks may be trimmed off during constant propagation.
    std::set<uint32_t> m_downstreamSpecMTaskIDs;

    ExecMTask* m_partnerSpecMTaskp
        = nullptr;  // Partner speculative MTask. Ie. if this MTask is the true branch of a
                    // speculation, this pointer will point to the 'false' branch mtask.
    VL_UNCOPYABLE(ExecMTask);

    string specSuffix() const {
        if (m_spec) {
            return (m_spec.value() ? "__SPEC__t_" : "__SPEC__f_") + cvtToStr(m_speculatesMtaskId);
        } else {
            return "";
        }
    }

public:
    ExecMTask(V3Graph* graphp, AstMTaskBody* bodyp, uint32_t id)
        : AbstractMTask{graphp}
        , m_bodyp{bodyp}
        , m_id{id} {}
    void speculative(bool specBranch, unsigned speculatesMtaskId) {
        m_spec = specBranch;
        m_speculatesMtaskId = speculatesMtaskId;
    }
    ExecMTask* partnerSpecMTask() { return m_partnerSpecMTaskp; }
    void partnerSpecMTask(ExecMTask* partnerSpecMTaskp) {
        m_partnerSpecMTaskp = partnerSpecMTaskp;
    }

    void addDownstreamSpeculativeTask(uint32_t mtaskid) {
        m_downstreamSpecMTaskIDs.insert(mtaskid);
    }
    const std::set<uint32_t>& downstreamSpeculativeMTasks() const {
        return m_downstreamSpecMTaskIDs;
    }

    void addUpstreamSpeculativeDepMTasks(uint32_t id) {
        m_upstreamSpeculativeDepMTasks.insert(id);
    }
    const std::set<uint32_t>& upstreamSpeculativeDepMTasks() const {
        return m_upstreamSpeculativeDepMTasks;
    }
    ExecMTask* partnerSpecMTaskp() const { return m_partnerSpecMTaskp; }
    std::optional<bool> speculative() const { return m_spec; }
    unsigned specMTaskId() const { return m_speculatesMtaskId; }
    AstMTaskBody* bodyp() const { return m_bodyp; }
    virtual uint32_t id() const override { return m_id; }
    uint32_t priority() const { return m_priority; }
    void priority(uint32_t pri) { m_priority = pri; }
    uint32_t downstreamCost() const { return m_downstreamCost; }
    void downstreamCost(uint32_t cost) { m_downstreamCost = cost; }
    uint32_t downstreamCostIncl() const { return m_downstreamCost + cost(); }
    virtual uint32_t cost() const override { return m_cost; }
    void cost(uint32_t cost) { m_cost = cost; }
    void thread(uint32_t thread) { m_thread = thread; }
    void startTime(uint32_t startTime) { m_startTime = startTime; }
    uint32_t startTime() const { return m_startTime; }
    uint32_t endTime() const { return m_endTime; }
    void endTime(uint32_t endTime) { m_endTime = endTime; }
    uint32_t thread() const { return m_thread; }
    void packNextp(const ExecMTask* nextp) { m_packNextp = nextp; }
    const ExecMTask* packNextp() const { return m_packNextp; }
    bool threadRoot() const { return m_threadRoot; }
    void threadRoot(bool threadRoot) { m_threadRoot = threadRoot; }
    string cFuncName() const {
        // If this MTask maps to a C function, this should be the name
        return string("__Vmtask") + "__" + cvtToStr(m_id) + specSuffix();
    }
    virtual string name() const override { return string("mt") + cvtToStr(id()) + specSuffix(); }
    void dump(std::ostream& str) const {
        str << name() << "." << cvtToHex(this);
        if (priority() || cost()) str << " [pr=" << priority() << " c=" << cvtToStr(cost()) << "]";
        if (thread() != 0xffffffff) str << " th=" << thread();
        if (threadRoot()) str << " [ROOT]";
        if (packNextp()) str << " nx=" << packNextp()->name();
    }
};
inline std::ostream& operator<<(std::ostream& os, const ExecMTask& rhs) {
    rhs.dump(os);
    return os;
}

#endif  // Guard
