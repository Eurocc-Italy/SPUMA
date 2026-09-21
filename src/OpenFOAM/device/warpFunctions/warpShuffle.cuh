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
    Implement shfl_sync functions and define specializations
    for selected vectorSpace types (vector,tensor)

SourceFiles
    warpShuffle.cuh

\*---------------------------------------------------------------------------*/

#ifndef Foam_device_warpShuffle_H
#define Foam_device_warpShuffle_H

#include "scalar.H"
#include "label.H"
#include "VectorSpace.H"

namespace Foam
{

namespace cuda
{

template <typename T>
__device__ __inline__ T __shfl_down_sync(unsigned mask, T val, uint offset, int width = warpSize)
{
    return ::__shfl_down_sync(mask, val, offset, width);
}

template <typename T>
__device__ __inline__ T __shfl_up_sync(unsigned mask, T val, uint offset, int width = warpSize)
{
    return ::__shfl_up_sync(mask, val, offset, width);
}
template <typename T>
__device__ __inline__ T __shfl_sync(unsigned mask, T val, int srcLane, int width = warpSize)
{
    return ::__shfl_sync(mask, val, srcLane, width);
}


#define SHFL_FUNCS_FOR_VECTORSPACE(Type)                                  \
__device__ __inline__ Type __shfl_down_sync(                              \
unsigned mask,                                                            \
Type &val,                                                                \
uint offset,                                                              \
int width = warpSize)                                                     \
{                                                                         \
Type tmp;                                                                 \
_Pragma("unroll") for (direction cmp = 0; cmp < Type::nComponents; ++cmp) \
{                                                                         \
    tmp[cmp] = ::__shfl_down_sync(mask, val[cmp], offset, width);         \
}                                                                         \
                                                                          \
return tmp;                                                               \
}                                                                         \
                                                                          \
__device__ __inline__ Type __shfl_up_sync(                                \
unsigned mask,                                                            \
Type &val,                                                                \
uint offset,                                                              \
int width = warpSize)                                                     \
{                                                                         \
Type tmp;                                                                 \
_Pragma("unroll") for (direction cmp = 0; cmp < Type::nComponents; ++cmp) \
{                                                                         \
    tmp[cmp] = ::__shfl_up_sync(mask, val[cmp], offset, width);           \
}                                                                         \
                                                                          \
return tmp;                                                               \
}                                                                         \
__device__ __inline__ Type __shfl_sync(                                   \
unsigned mask,                                                            \
Type &val,                                                                \
int srcLane,                                                              \
int width = warpSize)                                                     \
{                                                                         \
Type tmp;                                                                 \
_Pragma("unroll") for (direction cmp = 0; cmp < Type::nComponents; ++cmp) \
{                                                                         \
    tmp[cmp] = ::__shfl_sync(mask, val[cmp], srcLane, width);             \
}                                                                         \
return tmp;                                                               \
}

SHFL_FUNCS_FOR_VECTORSPACE(vector)
SHFL_FUNCS_FOR_VECTORSPACE(Vector<label>)
SHFL_FUNCS_FOR_VECTORSPACE(tensor)
SHFL_FUNCS_FOR_VECTORSPACE(Tensor<label>)

#undef SHFL_FUNCS_FOR_VECTORSPACE

}; // namespace cuda

}; // namespace Foam

#endif // Foam_device_warpShuffle_H
