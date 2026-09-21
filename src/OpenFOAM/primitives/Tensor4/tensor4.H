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

#ifndef Foam_tensor4_H
#define Foam_tensor4_H

#include "Tensor4.H"
#include "scalar.H"

namespace Foam
{

//! \class Foam::floatTensor4
//! \brief A Tensor4 of values with float precision
typedef Tensor4<float> floatTensor4;

//! \class Foam::doubleTensor4
//! \brief A Tensor4 of values with double precision
typedef Tensor4<double> doubleTensor4;

//! \class Foam::tensor4
//! \brief A Tensor4 with the configured scalar precision
typedef Tensor4<scalar> tensor4;

} // End namespace Foam

#endif

// ************************************************************************* //