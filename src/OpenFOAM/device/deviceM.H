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

Description
    Collection of device MACROs

SourceFiles
    deviceM.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_device_M_H
#define Foam_device_M_H

#ifdef have_cuda
    #include <cuda_runtime_api.h>
#elif have_hip
    #include <hip/hip_runtime.h>
#endif

#if defined(have_cuda) || defined(have_hip)
    #define FOAM_DEVICE __host__ __device__
    #define FOAM_DEVICE_I __host__ __device__ __forceinline__
    #define FOAM_DEVICE_V __device__
    // NVIDIA platform has a 32 bit mask (each bit represent a lane in the warp)
    // while AMD has a 64 bit mask for historical reason:
    // AMD devices have 64 thread per warp for gfx9 and 32 for gfx10 and above.
    #if defined(have_cuda)
        #define FULL_WARP_MASK_V 0xFFFFFFFF
        #define WARP_MASK_T unsigned int
    #elif defined(have_hip)
        #define FULL_WARP_MASK_V 0xFFFFFFFFFFFFFFFFULL
        #define WARP_MASK_T unsigned long long
    #endif
#else
    #define FOAM_DEVICE
    #define FOAM_DEVICE_I inline
    #define FOAM_DEVICE_V
    #define FULL_WARP_MASK_V
    #define WARP_MASK_T
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
