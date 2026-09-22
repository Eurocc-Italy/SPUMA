/*---------------------------------------------------------------------------*\
  *      .  *_______ * ______ .  __ *  __ * ___ .___    .  ___ .   *  .     *
    *  .    /       | |   _  \  |  |  |  | |   \/   | *   /   \ *   .    *   .
 *    .  * .\   (---*.|  |_)  |.|  |  |  |*|  \  /  |. * /  *  \  .  *     *
 =^^=^^==^^^=\   \^=^=|   ___/=^|  |^=|  |=|  |\/|  |^^=/  /=\  \^=^=^^===^^^=
 0  o  O  o---)   \ 0 |  |   0  |  o--o  |o|  |  |  | o/  _____  \ 0   o  O
     0    |_______/   |__| o   o \______/  |__| 0|__| /__/  o  \__\   o
  O   o  o        0  o      0   O        o    o       O  o     0   o    0  o
-------------------------------------------------------------------------------
    Copyright (C) 2025-2026 Cineca
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

#ifndef Foam_hip_executor_cu
#define Foam_hip_executor_cu

#include "hipExecutor.hpp"
#include "deviceInits.H"
#include "deviceM.H"
#include "hipError.hpp"
#include "deviceUtils.H"
#include "parallelForKernel.H"
#include "reductionKernel.H"
#include "recursivePrefixSumKernel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<typename F>
void Foam::hipExecutor::_backendFor(F& lambda, const label& size)
{
    if (size <= 0)
        return;

    const label nThreadsPerBlock = hipDeviceInit::getNumberOfThreadsPerBlock();
    const label numBlocks = device::setNumBlocks(size, nThreadsPerBlock);

    hipLaunchKernelGGL
    (
        Foam::device::lambdaKernel<F>,
        numBlocks,
        nThreadsPerBlock,
        0,
        0,
        lambda,
        size
    );

    hipDeviceSynchronize();
    CHECK_LAST_HIP_ERROR();
};

template<typename F>
void Foam::hipExecutor::_backendSerialFor(F& lambda, const label& size)
{
    if (size <= 0)
        return;

    hipLaunchKernelGGL
    (
        Foam::device::lambdaKernel<F>,
        1,
        1,
        0,
        0,
        lambda,
        size
    );

    hipDeviceSynchronize();
    CHECK_LAST_HIP_ERROR();
};

template <typename F,typename resultT>
void Foam::hipExecutor::_backendReductionSum(
    F& lambda,
    resultT* const __restrict__ result,
    const label& size
)
{
    if (size <= 0) return;

    resultT* dPtrResult = static_cast<resultT*>
    (
        Spuma::MemoryPool::getInstance()->allocate(sizeof(resultT))
    );

    Spuma::MemoryPool::getInstance()->memSet
    (
        (void*) dPtrResult,
        (const void*) result,
        sizeof(resultT),
        sizeof(resultT)
    );

    // create mutex
    Foam::Mutex mutex;

    const label nThreadsPerBlock = hipDeviceInit::getNumberOfThreadsPerBlock();
    const label numBlocks = device::setTreeReduceNumBlocks(size, nThreadsPerBlock);
    const label nStreamingMultiprocessors = hipDeviceInit::getNumberOfStreamingMultiprocessors();

    const int maxbytes = hipDeviceInit::getSharedMemoryPerBlock();

    // declare that this kernel can use up to maxbytes of dynamically allocated shared memory
    CHECK_HIP_ERROR
    (
        hipFuncSetAttribute
        (
            reinterpret_cast<const void*>(Foam::device::reductionLambdaSumKernel<resultT,F>),
            hipFuncAttributeMaxDynamicSharedMemorySize,
            maxbytes
        )
    );

    // NOTE: for multiple template arguments in the kernel use parentheses or the
    //       preprocessor cannot parse the hipLaunchKernelGGL macro correctly
    hipLaunchKernelGGL
    (
        (Foam::device::reductionLambdaSumKernel<resultT,F>),
        (numBlocks + nStreamingMultiprocessors - 1)/nStreamingMultiprocessors,
        nThreadsPerBlock,
        maxbytes,
        0,
        dPtrResult,
        lambda,
        mutex.getMutex(),
        size
    );

    hipDeviceSynchronize();
    CHECK_LAST_HIP_ERROR();

    Spuma::MemoryPool::getInstance()->copyOut(dPtrResult, result, sizeof(resultT));
    Spuma::MemoryPool::getInstance()->free(dPtrResult);
};

template <typename F, typename Op, typename resultT>
void Foam::hipExecutor::_backendReductionCompare
(
    F& lambda,
    Op& op,
    resultT* const __restrict__ result,
    const label& size
)
{
    if (size <= 0) return;

    resultT* dPtrResult = static_cast<resultT*>
    (
        Spuma::MemoryPool::getInstance()->allocate(sizeof(resultT))
    );

    Spuma::MemoryPool::getInstance()->memSet
    (
        (void*) dPtrResult,
        (const void*) result,
        sizeof(resultT),
        sizeof(resultT)
    );

    // create mutex
    Foam::Mutex mutex;

    const label nThreadsPerBlock = hipDeviceInit::getNumberOfThreadsPerBlock();
    const label numBlocks = device::setTreeReduceNumBlocks(size, nThreadsPerBlock);
    const label nStreamingMultiprocessors = hipDeviceInit::getNumberOfStreamingMultiprocessors();

    const label maxbytes = hipDeviceInit::getSharedMemoryPerBlock();

    // declare that this kernel can use up to maxbytes of dynamically allocated shared memory
    CHECK_HIP_ERROR
    (
        hipFuncSetAttribute
        (
            reinterpret_cast<const void*>(Foam::device::reductionLambdaCompareKernel<resultT,F,Op>),
            hipFuncAttributeMaxDynamicSharedMemorySize,
            maxbytes
        )
    );

    // NOTE: for multiple template arguments in the kernel use parentheses or the
    //       preprocessor cannot parse the hipLaunchKernelGGL macro correctly
    hipLaunchKernelGGL
    (
        (Foam::device::reductionLambdaCompareKernel<resultT,F,Op>),
        (numBlocks + nStreamingMultiprocessors - 1)/nStreamingMultiprocessors,
        nThreadsPerBlock,
        maxbytes,
        0,
        dPtrResult,
        lambda,
        op,
        mutex.getMutex(),
        size
    );

    hipDeviceSynchronize();
    CHECK_LAST_HIP_ERROR();

    Spuma::MemoryPool::getInstance()->copyOut(dPtrResult, result, sizeof(resultT));
    Spuma::MemoryPool::getInstance()->free(dPtrResult);
};

template<typename valueT>
void Foam::hipExecutor::_backendExclusiveSum
(
    valueT* const __restrict__ output,
    valueT* const __restrict__ input,
    const label& size
)
{
    if (size <= 0) return;

    constexpr size_t blockSize = 1024; // hardcode block size, best performance from tests
    constexpr uint element_per_thread = Foam::device::threadLocal<valueT>::element_per_thread;
    constexpr size_t element_per_block = element_per_thread * blockSize;

    const size_t num_blocks = (size + element_per_block - 1) / element_per_block;
    const size_t smem_size = element_per_block * sizeof(valueT);
    valueT* block_sum = nullptr;
    if(num_blocks > 1)
    {
        block_sum = static_cast<valueT*>
        (
            Spuma::MemoryPool::getInstance()->allocate(num_blocks*sizeof(valueT))
        );

        Spuma::MemoryPool::getInstance()->memSet
        (
            (void*) block_sum,
            0,
            num_blocks*sizeof(valueT)
        );
    }


    Foam::device::prefixSumKernel<valueT><<<num_blocks,blockSize,smem_size>>>
    (
        output,
        input,
        size,
        block_sum
    );

    hipDeviceSynchronize();
    CHECK_LAST_HIP_ERROR();

    // recursive block scan of the partial sums
    if(num_blocks > 1)
    {
        Foam::hipExecutor::_backendExclusiveSum<valueT>(block_sum,block_sum,num_blocks);

        Foam::device::addBlockKernel<valueT><<<num_blocks,blockSize>>>(
            output,
            block_sum,
            size,
            element_per_block
        );
        hipDeviceSynchronize();
        CHECK_LAST_HIP_ERROR();
        Spuma::MemoryPool::getInstance()->free(block_sum);
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
