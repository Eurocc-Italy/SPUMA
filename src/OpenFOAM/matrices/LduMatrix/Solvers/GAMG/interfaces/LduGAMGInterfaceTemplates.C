/*---------------------------------------------------------------------------*\\
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

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LduGAMGInterface.H"

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::LduGAMGInterface::interfaceInternalField
(
    const UList<Type>& internalData
) const
{
    auto tresult = tmp<Field<Type>>::New(size());
    interfaceInternalField(internalData, tresult.ref());
    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::LduGAMGInterface::interfaceInternalField
(
    const UList<Type>& internalData,
    const labelUList& faceCells
) const
{
    auto tresult = tmp<Field<Type>>::New(faceCells.size());
    auto& result = tresult.ref();

    forAll(result, elementi)
    {
        result[elementi] = internalData[faceCells[elementi]];
    }

    return tresult;
}


template<class Type>
void Foam::LduGAMGInterface::interfaceInternalField
(
    const UList<Type>& internalData,
    List<Type>& result
) const
{
    result.resize(size());

    forAll(result, elementi)
    {
        result[elementi] = internalData[faceCells_[elementi]];
    }
}


// ************************************************************************* //