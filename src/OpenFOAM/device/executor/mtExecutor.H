/*---------------------------------------------------------------------------*\
  *      .  *_______ * ______ .  __ *  __ * ___ .___    .  ___ .   *  .     *
    *  .    /       | |   _  \  |  |  |  | |   \/   | *   /   \ *   .    *   .
 *    .  * .\   (---*.|  |_)  |.|  |  |  |*|  \  /  |. * /  *  \  .  *     *
 =^^=^^==^^^=\   \^=^=|   ___/=^|  |^=|  |=|  |\/|  |^^=/  /=\  \^=^=^^===^^^=
 0  o  O  o---)   \ 0 |  |   0  |  o--o  |o|  |  |  | o/  _____  \ 0   o  O
     0    |_______/   |__| o   o \______/  |__| 0|__| /__/  o  \__\   o
  O   o  o        0  o      0   O        o    o       O  o     0   o    0  o
-------------------------------------------------------------------------------
    Copyright (C) 2026 Cineca
-------------------------------------------------------------------------------
License
    This file is part of SPUMA.

    SPUMA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SPUMA is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with SPUMA.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::mtExecutor

Description
    CPU multi-thread executor backend.
    Spawn local thread pool for each operation

SourceFiles
    mtExecutor.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_cpu_mt_executor_H
#define Foam_cpu_mt_executor_H


#include <thread>
#include <mutex>
#include <vector>
#include "executor.H"
#include "deviceInits.H"
#include "zero.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mtExecutor;

/*---------------------------------------------------------------------------*\
                          Class mtExecutor Declaration
\*---------------------------------------------------------------------------*/

class mtExecutor
:
    public executor<mtExecutor>
{

public:

    template<typename F>
    void _backendFor(F& lambda, const label& size)
    {
        if (size <= 0)
            return;

        const label nThreads = mtDeviceInit::getNumberOfThreadsPerBlock();
        std::vector<std::thread> workers;
        workers.reserve(nThreads);

        const int gridSize = (size + nThreads - 1) / nThreads;
        int grid_idx{0};
        for(label j = 0 ; j < nThreads; j++)
        {
            workers.emplace_back(
                std::thread([&,grid_idx](){
                for (label i = grid_idx; i < grid_idx + gridSize && i < size; i++)
                {
                    lambda(i);
                }
                }
                )
            );
            grid_idx+=gridSize;
        }

        // sync
        for (auto &&worker : workers)
        {
            worker.join();
        }
    }

    template<typename F>
    void _backendSerialFor(F& lambda, const label& size)
    {
        for (label i = 0; i < size; i++)
        {
            lambda(i);
        }
    }

    template<typename F, typename resultT>
    void _backendReductionSum
    (
        F& lambda,
        resultT* const __restrict__ result,
        const label& size
    )
    {
        // for now spawn local threads

        const label nThreads = mtDeviceInit::getNumberOfThreadsPerBlock();

        std::vector<std::thread> workers;
        workers.reserve(nThreads);

        const int gridSize = (size + nThreads - 1) / nThreads;
        int grid_idx{0};

        std::mutex reduction_mutex;

        for(label j = 0 ; j < nThreads; j++)
        {
            workers.emplace_back(
                std::thread([&,grid_idx](){
                    resultT local_result(Zero);

                    for (label i = grid_idx; i < grid_idx + gridSize && i<size; i++)
                    {
                        local_result += lambda(i);
                    }

                    {
                        std::lock_guard<std::mutex> lock(reduction_mutex);
                        *result += local_result;
                    }
                    }
                )
            );
            grid_idx+=gridSize;
        }
        
        // sync
        for (auto &&worker : workers)
        {
            worker.join();
        }
    }

    template<typename F,typename Op, typename resultT>
    void _backendReductionCompare
    (
        F& lambda,
        Op& op,
        resultT* const __restrict__ result,
        const label& size
    )
    {
        const label nThreads = mtDeviceInit::getNumberOfThreadsPerBlock();

        std::vector<std::thread> workers;
        workers.reserve(nThreads);

        const int gridSize = (size + nThreads - 1) / nThreads;
        int grid_idx{0};

        std::mutex reduction_mutex;

        for(label j = 0 ; j < nThreads; j++)
        {
            workers.emplace_back(
                std::thread([&,grid_idx](){
                    resultT local_result(*result);
                    for (label i = grid_idx; i < grid_idx + gridSize && i<size; i++)
                    {
                        local_result = op(lambda(i),local_result);
                    }

                    {
                        std::lock_guard<std::mutex> lock(reduction_mutex);
                        *result = op(local_result,*result);
                    }

                    }
                )
            );
            grid_idx+=gridSize;
        }
        
        // sync
        for (auto &&worker : workers)
        {
            worker.join();
        }
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
