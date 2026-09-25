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

Description
    Implement parallel gpu primitive kernels for:
    - reduction operations that apply a lambda function to each element and
        sum the results
    - reduction operations that apply a lambda function to elements and
        compute a binary comparison-based reductions operation

Reference:
    \verbatim
        Reduction algorithm reference:
            Harris, Mark. "Optimizing parallel reduction in
            CUDA." Nvidia developer technology 2.4 (2007): 70.
    \endverbatim

SourceFiles
    lambdaKernels.H

\*---------------------------------------------------------------------------*/
#ifndef Foam_device_reductionKernels_H
#define Foam_device_reductionKernels_H

#include "label.H"
#include "vector.H"
#include "scalar.H"
#include "tensor.H"
#include "deviceM.H"
#include "warpReduce.H"
#include "mutex.H"
#include "spinLock.H"
#include "sharedMemory.H"

// * * * * * * * * * * * * * * * * Kernels  * * * * * * * * * * * * * * //

namespace Foam
{

namespace device
{

template <typename resultType, typename F>
__global__
void reductionLambdaSumKernel
(
    resultType* const __restrict__ result,
    F lambda,
    int* const __restrict__ mutex,
    const label size
)
{
    unsigned int id  = blockIdx.x * (2*blockDim.x) + threadIdx.x;
    const unsigned int tid = threadIdx.x;
    const label gsize = size;
    const unsigned int blockSize = blockDim.x;
    const unsigned int gridSize = blockDim.x*2*gridDim.x;

    SharedMemory<resultType> smem;
    resultType *sdata = smem.getPointer();
    memset(&sdata[tid], 0, sizeof(resultType));

    __syncthreads();

    // grid-wise reduction step
    // load from gloabl memory + gsize/gridSize step of reduction
    while (id < gsize)
    {
        // handle loop_len not multiple of blockSize
        if (id + blockSize < gsize)
        {
            sdata[tid] += lambda(id) + lambda(id+blockSize);
        }
        else
        {
            sdata[tid] += lambda(id);
        }
        id += gridSize;
    }

    __syncthreads();

    // block wise reduction steps
    if (blockSize >= 1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
    if (blockSize >= 512)  { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >= 256)  { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >= 128)  { if (tid < 64)  { sdata[tid] += sdata[tid + 64];  } __syncthreads(); }

    // warp wise reduction step
    //in warp: syncthread guaranteed
    if (tid < 32)
    {
        warpReduceNoVolatile<resultType>(sdata, tid, blockSize);
    }

    __syncthreads();

    if constexpr(
        std::is_same_v<resultType,Foam::scalar>     ||
        std::is_same_v<resultType,Foam::label>      ||
        std::is_same_v<resultType,Foam::solveScalar>
    )
    {
        if(tid == 0)
        {
            atomicAdd(result,sdata[0]);
        }
    }
    else if constexpr(std::is_same_v<resultType,Foam::vector> || std::is_same_v<resultType,Foam::tensor>)
    {
        // parallel atomics for components of vector/tensor types
        if (tid < resultType::nComponents)
        {
            atomicAdd(reinterpret_cast<typename resultType::cmptType*>(result) + tid, sdata[0][tid]);
        }
    }
    else
    {
        // note for high number of block too much contention of the mutex!
        if(tid == 0)
        {
            spinLock::lock(mutex);
            __threadfence();
            *result += sdata[0];
            __threadfence();
             spinLock::unlock(mutex);
        }
    }
};

// room for improvement
template <typename resultType, typename F, typename Op>
__global__
void reductionLambdaCompareKernel
(
    resultType* const __restrict__ result,
    F lambda,
    Op op,
    int* const __restrict__ mutex,
    const label size
)
{
    unsigned int id  = blockIdx.x * (2*blockDim.x) + threadIdx.x;
    const unsigned int tid = threadIdx.x;
    const label gsize = size;
    const unsigned int blockSize = blockDim.x;
    const unsigned int gridSize = blockDim.x*2*gridDim.x;

    SharedMemory<resultType> smem;
    resultType *sdata = smem.getPointer();
    resultType tmp;
    sdata[tid] = *result;

    __syncthreads();

    // grid-wise reduction step
    // load from gloabl memory + gsize/gridSize step of reduction
    while (id < gsize)
    {
        // handle loop_len not multiple of blockSize
        if (id+blockSize < gsize)
        {
            tmp = op(lambda(id), lambda(id + blockSize));
        }
        else
        {
            tmp = lambda(id);
        }
        sdata[tid]=op(tmp,sdata[tid]);
        id += gridSize;
    }

    __syncthreads();

    // block wise reduction steps
    if (blockSize >= 1024)
    {
        if (tid < 512)
        {
            tmp = op(sdata[tid], sdata[tid + 512]);
            __threadfence_block();
            sdata[tid] = tmp;
        }
        __syncthreads();
    }
    if (blockSize >= 512)
    {
        if (tid < 256)
        {
            tmp = op(sdata[tid], sdata[tid + 256]);
            __threadfence_block();
            sdata[tid] = tmp;
        }
        __syncthreads();
    }
    if (blockSize >= 256)
    {
        if (tid < 128)
        {
            tmp = op(sdata[tid], sdata[tid + 128]);
            __threadfence_block();
            sdata[tid] = tmp;
        }
        __syncthreads();
    }
    if (blockSize >= 128)
    {
        if (tid < 64)
        {
            tmp = op(sdata[tid], sdata[tid + 64]);
            __threadfence_block();
                sdata[tid] = tmp;
        }
        __syncthreads();
    }

    // warp wise reduction step
    if (tid < 32)
    {
        warpReduceCompareNoVolatile<resultType,Op>(sdata,op,tid,blockSize);
    }

    __syncthreads();
    // note for high number of block too much contention of the mutex!
    if(tid == 0)
    {
        spinLock::lock(mutex);
        __threadfence();
        tmp = op(sdata[0],*result);
        __threadfence();
        *result = tmp;
        __threadfence();
        spinLock::unlock(mutex);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace device

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
