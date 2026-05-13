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

Struct
    Foam::mtAtomic

Description
    cpu multi thread Atomic backend

SourceFiles
    mtAtomic.H

\*---------------------------------------------------------------------------*/


#ifndef Foam_cpu_mt_Atomic_H
#define Foam_cpu_mt_Atomic_H

#include <atomic>
#include "Atomic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Struct mtAtomic Declaration
\*---------------------------------------------------------------------------*/


struct mtAtomic
: Atomic<mtAtomic>
{
    template<class T>
    struct atomicPlusEqOp
    {
        void operator()(T& x, const T& y) const
        {
            std::atomic_ref<T>(x) += y;
        }
    };

    template<class T>
    struct atomicMaxEqOp
    {
        void operator()(T& x, const T& y) const
        {
            _backendAtomicMax(x, y);
        }
    };

    template<class T>
    struct atomicMinEqOp
    {
        void operator()(T& x, const T& y) const
        {
            _backendAtomicMin(x, y);
        }
    };

    static void _backendAtomicAdd(scalar& x, const scalar& y)
    {
        std::atomic_ref<scalar>(x) += y;
    }

#ifdef WM_SPDP
    static void _backendAtomicAdd(solveScalar& x, const solveScalar& y)
    {
        std::atomic_ref<solveScalar>(x) += y;
    }
#endif

    static void _backendAtomicAdd(label& x, const label& y)
    {
        std::atomic_ref<label>(x) += y;
    }

    static void  _backendAtomicMax(scalar& x, const scalar& y)
    {
        std::atomic_ref<scalar> tmp(x);
        auto t = tmp.load();
        while(max(y, t) != t)
        {
            if (tmp.compare_exchange_weak(t,y)) break;
        }
    }

#ifdef WM_SPDP
    static void  _backendAtomicMax(solveScalar& x, const solveScalar& y)
    {
        std::atomic_ref<solveScalar> tmp(x);
        auto t = tmp.load();
        while(max(y, t) != t)
        {
            if (tmp.compare_exchange_weak(t,y)) break;
        }
    }
#endif

    static void  _backendAtomicMax(label& x, const label& y)
    {
        std::atomic_ref<label> tmp(x);
        auto t = tmp.load();
        while(max(y, t) != t)
        {
            if (tmp.compare_exchange_weak(t,y)) break;
        }
    }

    static void  _backendAtomicMin(scalar& x, const scalar& y)
    {
        std::atomic_ref<scalar> tmp(x);
        auto t = tmp.load();
        while(min(y, t) != t)
        {
            if (tmp.compare_exchange_weak(t,y)) break;
        }
    }

#ifdef WM_SPDP
    static void  _backendAtomicMin(solveScalar& x, const solveScalar& y)
    {
        std::atomic_ref<solveScalar> tmp(x);
        auto t = tmp.load();
        while(min(y, t) != t)
        {
            if (tmp.compare_exchange_weak(t,y)) break;
        }
    }
#endif

    static void  _backendAtomicMin(label& x, const label& y)
    {
        std::atomic_ref<label> tmp(x);
        auto t = tmp.load();
        while(min(y, t) != t)
        {
            if (tmp.compare_exchange_weak(t,y)) break;
        }
    }

    template<class Form, class Cmpt, direction Ncmpts>
    static void _backendAtomicAdd
    (
        VectorSpace<Form, Cmpt, Ncmpts>& vs1,
        const VectorSpace<Form, Cmpt, Ncmpts>& vs2
    )
    {
        VectorSpaceOps<Ncmpts,0>::eqOp(vs1, vs2, atomicPlusEqOp<Cmpt>());
    }

    template<class Form, class Cmpt, direction Ncmpts>
    static void _backendAtomicMax
    (
        VectorSpace<Form, Cmpt, Ncmpts>& vs1,
        const VectorSpace<Form, Cmpt, Ncmpts>& vs2
    )
    {
        VectorSpaceOps<Ncmpts,0>::eqOp(vs1, vs2, atomicMaxEqOp<Cmpt>());
    }

    template<class Form, class Cmpt, direction Ncmpts>
    static void _backendAtomicMin
    (
        VectorSpace<Form, Cmpt, Ncmpts>& vs1,
        const VectorSpace<Form, Cmpt, Ncmpts>& vs2
    )
    {
        VectorSpaceOps<Ncmpts,0>::eqOp(vs1, vs2, atomicMinEqOp<Cmpt>());
    }

    static label _backendAtomicCAS
    (
        label& x,
        const label& compare,
        const label& y
    )
    {
        return std::atomic_ref<label>(x).compare_exchange_weak(const_cast<label&>(compare),y);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
