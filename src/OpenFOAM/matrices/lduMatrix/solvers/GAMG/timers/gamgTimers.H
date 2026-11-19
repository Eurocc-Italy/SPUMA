/*---------------------------------------------------------------------------*\
    Added (parallelsmoothers): scoped timers for GAMGSolver, the ones of our
    GAMGTimed (parallelsmoothers gamg/src/gamgVecTimers.H, same region names),
    on the region profiler in sdv_tracing_local.H (adapted from
    parallelsmoothers smoothers/tracing). Every region ends up in
    profiling_results_rank<r>.csv as "gv_<name>" (whole run, one row) or
    "gv_<name>@n<nCells>" (one row per multigrid level). Regions may nest;
    each one counts its own wall time.
    On by default; "-gamg_profiling no" on the command line turns every
    region into a no-op: nothing is registered, timed or written.
\*---------------------------------------------------------------------------*/

#ifndef gamgTimers_H
#define gamgTimers_H

#include "sdv_tracing_local.H"
#include "gamgProfiling.H"
#include "UPstream.H"

#include <map>
#include <utility>

namespace Foam
{
namespace gv
{

class Timer
{
    int id_;

public:

    explicit Timer(int id) : id_(id) { trace_begin_id(id_); }
    ~Timer() { stop(); }
    void stop() { if (id_ >= 0) { trace_end_id(id_); id_ = -1; } }
    Timer(const Timer&) = delete;
    Timer& operator=(const Timer&) = delete;
};

//- Region id of "name"; the first call records the MPI rank, which names
//  the CSV file written at exit (profiling_results_rank<r>.csv)
inline int regionId(const char* name)
{
    if (trace_rank == -2)
    {
        trace_rank = UPstream::parRun() ? UPstream::myProcNo() : -1;
    }
    return trace_region_id(name);
}

//- Per-level region id, looked up once per (name, nCells)
inline int levelId(const char* name, long nCells)
{
    static std::map<std::pair<const char*, long>, int> cache;
    const auto key = std::make_pair(name, nCells);
    const auto it = cache.find(key);
    if (it != cache.end()) return it->second;
    char full[REGION_NAME_LEN];
    snprintf(full, sizeof full, "%s@n%ld", name, nCells);
    const int id = regionId(full);
    cache[key] = id;
    return id;
}

} // End namespace gv
} // End namespace Foam

//- Whole-run region "gv_<name>"; lives to the end of the enclosing scope
#define GV_TIMER(var, name) \
    static const int var##_id_ = \
        Foam::gv::enabled ? Foam::gv::regionId("gv_" name) : -1; \
    Foam::gv::Timer var(var##_id_)

//- Per-level region "gv_<name>@n<nCells>"; lives to the end of the scope
#define GV_TIMER_LEVEL(var, name, nCells) \
    Foam::gv::Timer var \
    ( \
        Foam::gv::enabled ? Foam::gv::levelId("gv_" name, nCells) : -1 \
    )

#endif
