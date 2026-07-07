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

#ifndef Foam_LduGAMGInterfaceField_H
#define Foam_LduGAMGInterfaceField_H

#include "CyclicLduGAMGInterfaceField.H"
#include "ProcessorCyclicLduGAMGInterfaceField.H"
#include "ProcessorLduGAMGInterfaceField.H"

namespace Foam
{

template<class Type>
autoPtr<LduInterfaceField<Type>> newLduGAMGInterfaceField
(
    const LduGAMGInterface& coarseInterface,
    const LduInterfaceField<Type>& fineInterface
)
{
    if (isA<processorCyclicLduGAMGInterface>(coarseInterface))
    {
        return autoPtr<LduInterfaceField<Type>>
        (
            new ProcessorCyclicLduGAMGInterfaceField<Type>
            (
                refCast<const processorCyclicLduGAMGInterface>
                (
                    coarseInterface
                ),
                fineInterface
            )
        );
    }

    if (isA<processorLduGAMGInterface>(coarseInterface))
    {
        return autoPtr<LduInterfaceField<Type>>
        (
            new ProcessorLduGAMGInterfaceField<Type>
            (
                refCast<const processorLduGAMGInterface>(coarseInterface),
                fineInterface
            )
        );
    }

    if (isA<cyclicLduGAMGInterface>(coarseInterface))
    {
        return autoPtr<LduInterfaceField<Type>>
        (
            new CyclicLduGAMGInterfaceField<Type>
            (
                refCast<const cyclicLduGAMGInterface>(coarseInterface),
                fineInterface
            )
        );
    }

    FatalErrorInFunction
        << "Unsupported LduMatrix GAMG interface field combination. "
        << "Coarse interface type: " << coarseInterface.type()
        << ", fine interface field type: "
        << fineInterface.interfaceFieldType()
        << exit(FatalError);

    return nullptr;
}

} // End namespace Foam

#endif

// ************************************************************************* //