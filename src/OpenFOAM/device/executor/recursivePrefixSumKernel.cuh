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
    recursivePrefixSumKernel.cuh

\*---------------------------------------------------------------------------*/
#ifndef Foam_recursive_prefix_sum_kernel_H
#define Foam_recursive_prefix_sum_kernel_H

#include "warpShuffle.cuh"
#include "deviceM.H"
#include "label.H"


namespace Foam
{

namespace cuda
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template <typename T>
struct type_tuning
{
    typedef T value_type;
    typedef T reduce_type;

    static constexpr int element_per_thread = 1;
};

#define TYPE_TUNING_SPECIALIZATION(scalar_type, N)                                  \
template <>                                                                         \
struct type_tuning<scalar_type>                                                     \
{                                                                                   \
typedef scalar_type##N value_type;                                                  \
typedef scalar_type reduce_type;                                                    \
static_assert(sizeof(value_type) % sizeof(reduce_type) == 0,                        \
              "value_type must be a multiple of the size of reduce_type");          \
static constexpr int element_per_thread = sizeof(value_type) / sizeof(reduce_type); \
};

TYPE_TUNING_SPECIALIZATION(int,4)
TYPE_TUNING_SPECIALIZATION(float,4)
TYPE_TUNING_SPECIALIZATION(double,2)


template <typename T>
struct threadLocal
{    
    using value_type = typename type_tuning<T>::value_type;
    using reduce_type = typename type_tuning<T>::reduce_type;
    static constexpr int element_per_thread = type_tuning<T>::element_per_thread;

    value_type _value;
    reduce_type _reduce_value;

    explicit threadLocal() : _value{}, _reduce_value{} {};


    __device__ __forceinline__ void padded_load
    (
        size_t global_id,
        int size,
        T* const __restrict__ input,
        const T default_v
    )
    {
        _value = (global_id < size) ? input[global_id] : default_v;
    }

    __device__ __forceinline__ void padded_load
    (
        size_t global_id,
        int size,
        T* const __restrict__ input
    )
    {
        const T default_v{};
        _value = (global_id < size) ? input[global_id] : default_v;
    }
    
    __device__ __forceinline__ T reduce()
    {
        _reduce_value = _value;
        return _value;
    }

    __device__ __forceinline__ void exclusive_sum_store
    (
        size_t global_id,
        int size,
        T* const __restrict__ output,
        const T& value
    )
    {
        if (global_id < size) output[global_id] = value;
    }
};

#define THREAD_LOCAL_SPECIALIZATION(scalar_type)                                                                  \
template <>                                                                                                       \
struct threadLocal<scalar_type>                                                                                   \
{                                                                                                                 \
using value_type = typename type_tuning<scalar_type>::value_type;                                                 \
using reduce_type = typename type_tuning<scalar_type>::reduce_type;                                               \
static constexpr int element_per_thread = type_tuning<scalar_type>::element_per_thread;                           \
value_type _value;                                                                                                \
reduce_type _reduce_value;                                                                                        \
                                                                                                                  \
explicit threadLocal() : _value{}, _reduce_value{} {};                                                            \
                                                                                                                  \
__device__ __forceinline__ void padded_load(                                                                      \
    size_t global_id,                                                                                             \
    int size,                                                                                                     \
    reduce_type *const __restrict__ input,                                                                        \
    const reduce_type default_v)                                                                                  \
{                                                                                                                 \
    /* initialize _value with default_v for all elements*/                                                        \
    _Pragma("unroll") for (int i = 0; i < element_per_thread; ++i)                                                \
    {                                                                                                             \
        reinterpret_cast<reduce_type *>(&(_value))[i] = default_v;                                                \
    };                                                                                                            \
    const size_t thread_global_id = global_id * element_per_thread;                                               \
    if (thread_global_id + element_per_thread - 1 < size)                                                         \
    {                                                                                                             \
        _value = reinterpret_cast<const value_type *>(input)[thread_global_id / element_per_thread];              \
    }                                                                                                             \
    else                                                                                                          \
    {                                                                                                             \
        _Pragma("unroll") for (int i = 0; i < element_per_thread; ++i)                                            \
        {                                                                                                         \
            if (thread_global_id + i < size)                                                                      \
                reinterpret_cast<reduce_type *>(&(_value))[i] = input[thread_global_id + i];                      \
        }                                                                                                         \
    }                                                                                                             \
};                                                                                                                \
__device__ __forceinline__ void padded_load(size_t global_id, int size, reduce_type *const __restrict__ input)    \
{                                                                                                                 \
    _value = value_type{}; /*initialize _value with zeros for all elements*/                                      \
    const size_t thread_global_id = global_id * element_per_thread;                                               \
    if (thread_global_id + element_per_thread - 1 < size)                                                         \
    {                                                                                                             \
        _value = reinterpret_cast<const value_type *>(input)[thread_global_id / element_per_thread];              \
    }                                                                                                             \
    else                                                                                                          \
    {                                                                                                             \
        _Pragma("unroll") for (int i = 0; i < element_per_thread; ++i)                                            \
        {                                                                                                         \
            if (thread_global_id + i < size)                                                                      \
                reinterpret_cast<reduce_type *>(&(_value))[i] = input[thread_global_id + i];                      \
        }                                                                                                         \
    }                                                                                                             \
};                                                                                                                \
__device__ __forceinline__ reduce_type reduce()                                                                   \
{                                                                                                                 \
    _Pragma("unroll") for (int i = 0; i < element_per_thread; ++i)                                                \
    {                                                                                                             \
        _reduce_value += reinterpret_cast<reduce_type *>(&(_value))[i];                                           \
    }                                                                                                             \
    return _reduce_value;                                                                                         \
};                                                                                                                \
__device__ __forceinline__ void exclusive_sum_store(                                                              \
    size_t global_id,                                                                                             \
    int size,                                                                                                     \
    reduce_type *const __restrict__ output,                                                                       \
    const reduce_type value)                                                                                      \
{                                                                                                                 \
    value_type out_val{};                                                                                         \
    reinterpret_cast<reduce_type *>(&out_val)[0] = value;                                                         \
    _Pragma("unroll") for (int i = 1; i < element_per_thread; ++i)                                                \
    {                                                                                                             \
        reinterpret_cast<reduce_type *>(&out_val)[i] =                                                            \
            reinterpret_cast<reduce_type *>(&out_val)[i - 1] + reinterpret_cast<reduce_type *>(&(_value))[i - 1]; \
    };                                                                                                            \
    const size_t thread_global_id = global_id * element_per_thread;                                               \
    if (thread_global_id + element_per_thread - 1 < size)                                                         \
    {                                                                                                             \
        reinterpret_cast<value_type *>(output)[thread_global_id / element_per_thread] = out_val;                  \
    }                                                                                                             \
    else                                                                                                          \
    {                                                                                                             \
        _Pragma("unroll") for (int i = 0; i < element_per_thread; ++i)                                            \
        {                                                                                                         \
            if (thread_global_id + i < size)                                                                      \
                output[thread_global_id + i] = reinterpret_cast<reduce_type *>(&out_val)[i];                      \
        };                                                                                                        \
    }                                                                                                             \
}                                                                                                                 \
};

THREAD_LOCAL_SPECIALIZATION(int)
THREAD_LOCAL_SPECIALIZATION(float)
THREAD_LOCAL_SPECIALIZATION(double)

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
    const int element_per_block = blockSize * threadLocal<T>::element_per_thread;

    const T zero_value{};
    // load data
    threadLocal<T> thread_data;
    thread_data.padded_load(global_id, size, input, zero_value);
    const T init_value = thread_data.reduce();
    T value = init_value;

    // warp shfl scan upsweep
    // inter warp scan
    constexpr unsigned mask = FULL_WARP_MASK_V;
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

} // namespace cuda


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// host side function

template<typename valueT>
void Foam::cudaExecutor::_backendExclusiveSum
(
    valueT* const __restrict__ output,
    valueT* const __restrict__ input,
    const label& size
)
{
    if (size <= 0) return;

    constexpr size_t blockSize = 1024; // hardcode block size, best performance from tests
    constexpr uint element_per_thread = Foam::cuda::threadLocal<valueT>::element_per_thread;
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


    Foam::cuda::prefixSumKernel<valueT><<<num_blocks,blockSize,smem_size>>>
    (
        output,
        input,
        size,
        block_sum
    );

    cudaDeviceSynchronize();
    CHECK_LAST_CUDA_ERROR();
    
    // recursive block scan of the partial sums
    if(num_blocks > 1)
    {   
        Foam::cudaExecutor::_backendExclusiveSum<valueT>(block_sum,block_sum,num_blocks);
  
        Foam::cuda::addBlockKernel<valueT><<<num_blocks,blockSize>>>(
            output,
            block_sum,
            size,
            element_per_block
        );
        cudaDeviceSynchronize();
        CHECK_LAST_CUDA_ERROR();
        Spuma::MemoryPool::getInstance()->free(block_sum);
    }

};


} // namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
