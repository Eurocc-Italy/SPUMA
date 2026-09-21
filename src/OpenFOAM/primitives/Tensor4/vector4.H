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

#ifndef Foam_vector4_H
#define Foam_vector4_H

#include "Vector4.H"
#include "scalar.H"

namespace Foam
{

//! \class Foam::floatVector4
//! \brief A Vector4 of values with float precision
typedef Vector4<float> floatVector4;

//! \class Foam::doubleVector4
//! \brief A Vector4 of values with double precision
typedef Vector4<double> doubleVector4;

//! \class Foam::vector4
//! \brief A Vector4 with the configured scalar precision
typedef Vector4<scalar> vector4;

} // End namespace Foam

#endif

// ************************************************************************* //