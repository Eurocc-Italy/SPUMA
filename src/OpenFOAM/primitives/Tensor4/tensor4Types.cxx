/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Cineca
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

\*---------------------------------------------------------------------------*/

#include "tensor4.H"
#include "vector4.H"

#if defined(WM_DP)
template<>
const char* const Foam::Tensor4<float>::vsType::typeName = "floatTensor4";
template<>
const char* const Foam::Tensor4<double>::vsType::typeName = "tensor4";
template<>
const char* const Foam::Vector4<float>::vsType::typeName = "floatVector4";
template<>
const char* const Foam::Vector4<double>::vsType::typeName = "vector4";
#else
template<>
const char* const Foam::Tensor4<float>::vsType::typeName = "tensor4";
template<>
const char* const Foam::Tensor4<double>::vsType::typeName = "doubleTensor4";
template<>
const char* const Foam::Vector4<float>::vsType::typeName = "vector4";
template<>
const char* const Foam::Vector4<double>::vsType::typeName = "doubleVector4";
#endif

#define defineVector4Traits(Type, Prefix)                                    \
    template<>                                                               \
    const char* const Foam::Vector4<Type>::vsType::componentNames[] =         \
        {"x", "y", "z", "k"};                                             \
    template<>                                                               \
    const Foam::Vector4<Type> Foam::Vector4<Type>::vsType::zero               \
        (Foam::Vector4<Type>::uniform(0));                                    \
    template<>                                                               \
    const Foam::Vector4<Type> Foam::Vector4<Type>::vsType::one                \
        (Foam::Vector4<Type>::uniform(1));                                    \
    template<>                                                               \
    const Foam::Vector4<Type> Foam::Vector4<Type>::vsType::max                \
        (Foam::Vector4<Type>::uniform(Foam::Prefix##VGREAT));                 \
    template<>                                                               \
    const Foam::Vector4<Type> Foam::Vector4<Type>::vsType::min                \
        (Foam::Vector4<Type>::uniform(-Foam::Prefix##VGREAT));                \
    template<>                                                               \
    const Foam::Vector4<Type> Foam::Vector4<Type>::vsType::rootMax            \
        (Foam::Vector4<Type>::uniform(Foam::Prefix##ROOTVGREAT));             \
    template<>                                                               \
    const Foam::Vector4<Type> Foam::Vector4<Type>::vsType::rootMin            \
        (Foam::Vector4<Type>::uniform(-Foam::Prefix##ROOTVGREAT));

#define defineTensor4Traits(Type, Prefix)                                    \
    template<>                                                               \
    const char* const Foam::Tensor4<Type>::vsType::componentNames[] =         \
    {                                                                        \
        "xx", "xy", "xz", "xk", "yx", "yy", "yz", "yk",            \
        "zx", "zy", "zz", "zk", "kx", "ky", "kz", "kk"             \
    };                                                                       \
    template<>                                                               \
    const Foam::Tensor4<Type> Foam::Tensor4<Type>::vsType::zero               \
        (Foam::Tensor4<Type>::uniform(0));                                    \
    template<>                                                               \
    const Foam::Tensor4<Type> Foam::Tensor4<Type>::vsType::one                \
        (Foam::Tensor4<Type>::uniform(1));                                    \
    template<>                                                               \
    const Foam::Tensor4<Type> Foam::Tensor4<Type>::vsType::max                \
        (Foam::Tensor4<Type>::uniform(Foam::Prefix##VGREAT));                 \
    template<>                                                               \
    const Foam::Tensor4<Type> Foam::Tensor4<Type>::vsType::min                \
        (Foam::Tensor4<Type>::uniform(-Foam::Prefix##VGREAT));                \
    template<>                                                               \
    const Foam::Tensor4<Type> Foam::Tensor4<Type>::vsType::rootMax            \
        (Foam::Tensor4<Type>::uniform(Foam::Prefix##ROOTVGREAT));             \
    template<>                                                               \
    const Foam::Tensor4<Type> Foam::Tensor4<Type>::vsType::rootMin            \
        (Foam::Tensor4<Type>::uniform(-Foam::Prefix##ROOTVGREAT));

defineVector4Traits(float, floatScalar)
defineVector4Traits(double, doubleScalar)
defineTensor4Traits(float, floatScalar)
defineTensor4Traits(double, doubleScalar)

template<>
const Foam::Tensor4<Foam::label> Foam::Tensor4<Foam::label>::vsType::zero
(
    Foam::Tensor4<Foam::label>::uniform(0)
);

template<>
const Foam::Tensor4<Foam::label> Foam::Tensor4<Foam::label>::vsType::one
(
    Foam::Tensor4<Foam::label>::uniform(1)
);

template<>
const Foam::Vector4<Foam::label> Foam::Vector4<Foam::label>::vsType::zero
(
    Foam::Vector4<Foam::label>::uniform(0)
);

template<>
const Foam::Vector4<Foam::label> Foam::Vector4<Foam::label>::vsType::one
(
    Foam::Vector4<Foam::label>::uniform(1)
);

#undef defineVector4Traits
#undef defineTensor4Traits

// ************************************************************************* //