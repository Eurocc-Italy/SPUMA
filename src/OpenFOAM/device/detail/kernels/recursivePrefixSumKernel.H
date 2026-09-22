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
    Implement parallel exlusive prefix sum for primitive and
    general vectorspace data type.
    The partial block prefix sum is based on Blelloch's algorithm and uses
    a divide and conquer approach for array sizes greater than the blocksize.
    With respect to the reference cited there are different design choises
    implemented:
    - Increase ILP for primitive data types (int,float,double) via per thread
    partial prefix sum
    - Use shfl_sync functions to reduce shared memory footprint,
    synchrnonization points and memory latency via warp level prefix sum steps


Reference:
    \verbatim
        Base algorithm reference:
            Harris, Mark, Sengupta, S., & Owens, J. D. (2007).
            Parallel prefix sum (scan) with CUDA. GPU gems, 3(39), 851-876.
    \endverbatim

SourceFiles
    recursivePrefixSumKernel.H

\*---------------------------------------------------------------------------*/
#ifndef Foam_recursive_prefix_sum_kernel_H
#define Foam_recursive_prefix_sum_kernel_H

#include "threadLocal.H"
#include "warpShuffle.H"
#include "deviceM.H"
#include "label.H"


namespace Foam
{

namespace device
{

// * * * * * * * * * * * * * * * * Kernels  * * * * * * * * * * * * * * //


template <typename T>
__global__
void addBlockKernel(
    T* const __restrict__ output,
    T* const __restrict__ block_sum,
    const int size,
    const int element_per_block
)
{
    unsigned int id  = blockIdx.x *blockDim.x + threadIdx.x;
    const unsigned int gridSize = blockDim.x*gridDim.x;

    while(id < size)
    {
        const unsigned int logical_bid = id / (element_per_block);
        output[id] += block_sum[logical_bid];
        id += gridSize;
    }
};

template <typename T>
__global__
void __launch_bounds__(1024) prefixSumKernel(
    T* const __restrict__ output,
    T* const __restrict__ input,
    int size,
    T* const __restrict__ block_sum = nullptr
)
{
    extern __shared__ char smem_raw[];
    T* s_data = reinterpret_cast<T*>(smem_raw);

    const int tid = threadIdx.x;
    const int blockSize = blockDim.x;
    const size_t global_id = blockIdx.x * blockSize + tid;
    const int lane_id = tid % warpSize;
    const int warp_id = tid / warpSize;
    const int n_warps = blockSize / warpSize; // assumes blockSize is a multiple of warpSize
    const int bid = blockIdx.x;
   //const int element_per_block = blockSize * threadLocal<T>::element_per_thread;

    const T zero_value{};
    // load data
    threadLocal<T> thread_data;
    thread_data.padded_load(global_id, size, input, zero_value);
    const T init_value = thread_data.reduce();
    T value = init_value;


    // NVIDIA platform has a 32 bit mask (each bit represent a lane in the warp)
    // while AMD has a 64 bit mask for historical reason:
    // AMD devices have 64 thread per warp for gfx9 and 32 for gfx10 and above.
    #if defined(have_cuda)
        using mask_t = unsigned int;
        constexpr mask_t mask = 0xFFFFFFFF;
    #elif defined(have_hip)
        using mask_t = unsigned long long;
        constexpr mask_t mask = 0xFFFFFFFFFFFFFFFFULL;
    #endif

    // warp shfl scan upsweep
    // inter warp scan
    #pragma unroll
    for(int offset = 1; offset < warpSize ; offset *= 2)
    {
        T tmp = __shfl_up_sync(mask, value, offset);
        if (lane_id >= offset) {
            value += tmp;
        }
    }

    // store warp-wide partial sums
    if (lane_id == warpSize -1)
    {
        s_data[warp_id] = value;
    }
    __syncthreads();

    // first warp does the inter-warp pass
    if (warp_id == 0)
    {
        T warp_sum = (tid < n_warps) ? s_data[tid] : zero_value;
        const T prev_value = warp_sum;
        #pragma unroll
        for (int offset = 1; offset < warpSize; offset *= 2) {
            const T tmp = __shfl_up_sync(mask, warp_sum, offset);
            if (lane_id >= offset) {
                warp_sum += tmp;
            }
        }

        if (tid < n_warps) {
            s_data[tid] = warp_sum - prev_value;
        }

    }
    __syncthreads();

    // scan down sweep
    // add the inter-warp scan result to each warp
    value += s_data[warp_id];

    // block aggregate
    if (block_sum && tid == blockSize - 1) {
        // value is inclusive scan within this warp
        // s_data[warp_id] is exclusive prefix of previous warps
        block_sum[bid] = value;
    }

    // convert to exclusive
    value -= init_value;

    // store data
    thread_data.exclusive_sum_store(global_id, size, output, value);
};

} // namespace device


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
