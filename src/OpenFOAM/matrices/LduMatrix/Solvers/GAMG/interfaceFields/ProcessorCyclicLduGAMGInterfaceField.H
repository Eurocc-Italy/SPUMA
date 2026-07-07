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

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef Foam_ProcessorCyclicLduGAMGInterfaceField_H
#define Foam_ProcessorCyclicLduGAMGInterfaceField_H

#include "ProcessorLduGAMGInterfaceField.H"
#include "processorCyclicLduGAMGInterface.H"

namespace Foam
{

template<class Type>
class ProcessorCyclicLduGAMGInterfaceField
:
    public ProcessorLduGAMGInterfaceField<Type>
{
public:

    ProcessorCyclicLduGAMGInterfaceField
    (
        const processorCyclicLduGAMGInterface& coarseInterface,
        const LduInterfaceField<Type>& fineInterface
    )
    :
        ProcessorLduGAMGInterfaceField<Type>(coarseInterface, fineInterface)
    {}

    //- Construct with explicit transform settings (used by proc agglomeration)
    ProcessorCyclicLduGAMGInterfaceField
    (
        const processorCyclicLduGAMGInterface& coarseInterface,
        const bool doTransform,
        const int rank
    )
    :
        ProcessorLduGAMGInterfaceField<Type>
        (
            coarseInterface,
            doTransform,
            rank
        )
    {}
};

} // End namespace Foam

#endif

// ************************************************************************* //