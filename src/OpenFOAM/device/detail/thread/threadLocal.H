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
    Define a thread-local  structure for handling per-thread data
    with vector type specialization for selected data types. This is
    done to increase overall throughput by hiding memory latency and
    increase ILP (Instruction Level Parallelism).


SourceFiles
    threadLocal.H

\*---------------------------------------------------------------------------*/
#ifndef Foam_threadLocal_H
#define Foam_threadLocal_H

namespace Foam
{

namespace device
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

    explicit __device__ threadLocal() : _value{}, _reduce_value{} {};


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
explicit __device__ threadLocal() : _value{}, _reduce_value{} {};                                                 \
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


#undef THREAD_LOCAL_SPECIALIZATION
#undef TYPE_TUNING_SPECIALIZATION

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace device

} // namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

