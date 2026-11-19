/*---------------------------------------------------------------------------*\
    Added: on/off flag of the GAMG region timers (gamgTimers.H).
    On by default; argList sets it from the command-line option
    -gamg_profiling <bool> (e.g. -gamg_profiling no) before any solver runs.
\*---------------------------------------------------------------------------*/

#ifndef gamgProfiling_H
#define gamgProfiling_H

namespace Foam
{
namespace gv
{

//- True: GAMG region timers record and profiling_results*.csv is written
extern bool enabled;

} // End namespace gv
} // End namespace Foam

#endif
