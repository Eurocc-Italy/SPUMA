/*---------------------------------------------------------------------------*\
  *      .  *_______ * ______ .  __ *  __ * ___ .___    .  ___ .   *  .     *
    *  .    /       | |   _  \  |  |  |  | |   \/   | *   /   \ *   .    *   .
 *    .  * .\   (---*.|  |_)  |.|  |  |  |*|  \  /  |. * /  *  \  .  *     *
 =^^=^^==^^^=\   \^=^=|   ___/=^|  |^=|  |=|  |\/|  |^^=/  /=\  \^=^=^^===^^^=
 0  o  O  o---)   \ 0 |  |   0  |  o--o  |o|  |  |  | o/  _____  \ 0   o  O
     0    |_______/   |__| o   o \______/  |__| 0|__| /__/  o  \__\   o
  O   o  o        0  o      0   O        o    o       O  o     0   o    0  o
-------------------------------------------------------------------------------
    Copyright (C) 2025 Cineca
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

\*---------------------------------------------------------------------------*/

#include "MemoryPool.H"
#include "error.H"
#ifdef have_umpire
    #include "umpireMemoryPool.H"
#endif
#include "fixedSizeMemoryPool.H"
#include "dummyMemoryPool.H"
#include "memoryExecutors.H"


// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MemoryPool,  0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MemoryPool::MemoryPool():
    size_(0),
    allocatedSize_(0),
    unallocatedSize_(0),
    maxOccupancy_(0)
{};

Foam::MemoryPool::MemoryPool(const dictionary& dict):
    size_(0),
    allocatedSize_(0),
    unallocatedSize_(0),
    maxOccupancy_(0)
{};


// * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //

Foam::MemoryPool::~MemoryPool()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::MemoryPool* Foam::MemoryPool::New
(
    const word arch,
    const word type,
    const uint64_t size
)
{
    MemoryPool* instance = nullptr;

    if (type == "fixedSizeMemoryPool")
    {
        if (arch == "cpu")
        {
            instance = new fixedSizeMemoryPool<cpuMemoryExecutor>(size);
        }
#ifdef have_cuda
        else if (arch == "cuda")
        {
            instance = new fixedSizeMemoryPool<cudaMemoryExecutor>(size);
        }   
#elif have_hip
        else if (arch == "hip")
        {
            instance = new fixedSizeMemoryPool<hipMemoryExecutor>(size);
        }
#endif
        else
        {
            FatalErrorInFunction << arch << " architecture is not supported."
                << exit(FatalError);
        }
    }
    else if (type == "dummyMemoryPool")
    {
        if (arch == "cpu")
        {
            instance = new dummyMemoryPool<cpuMemoryExecutor>(size);
        }
#ifdef have_cuda
        else if (arch == "cuda")
        {
            instance = new dummyMemoryPool<cudaMemoryExecutor>(size);
        }
#elif have_hip
        else if (arch == "hip")
        {
            instance = new dummyMemoryPool<hipMemoryExecutor>(size);
        }
#endif
        else
        {
            FatalErrorInFunction << arch << " architecture is not supported."
                << exit(FatalError);
        }
    }
#ifdef have_umpire
    else if (type == "umpireMemoryPool")
    {
        if (arch == "cpu")
        {
            instance = new umpireMemoryPool<cpuMemoryExecutor>(size);
        }
#ifdef have_cuda
        else if (arch == "cuda")
        {
            instance = new umpireMemoryPool<cudaMemoryExecutor>(size);
        }
#elif have_hip
        else if (arch == "hip")
        {
            instance = new umpireMemoryPool<hipMemoryExecutor>(size);
        }
#endif
        else
        {
            FatalErrorInFunction << arch << " architecture is not supported."
                << exit(FatalError);
        }
    }
#endif
    else
    {
        FatalErrorInFunction
            << type << " does not exist. "
            << "Please use a different memory pool." << nl
            << abort(FatalError);
    }
    return instance;
}


// ************************************************************************* //
