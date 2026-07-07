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

Class
    Foam::processorCyclicLduGAMGInterface

Description
    Native GAMG agglomerated processor-cyclic interface for templated
    LduMatrix GAMG.

SourceFiles
    processorCyclicLduGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_processorCyclicLduGAMGInterface_H
#define Foam_processorCyclicLduGAMGInterface_H

#include "processorLduGAMGInterface.H"

namespace Foam
{

class processorCyclicLduGAMGInterface
:
    public processorLduGAMGInterface
{
    // Private Member Functions

        processorCyclicLduGAMGInterface
        (
            const processorCyclicLduGAMGInterface&
        ) = delete;
        void operator=(const processorCyclicLduGAMGInterface&) = delete;


public:

    TypeName("processorCyclic");


    // Constructors

        processorCyclicLduGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& localRestrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );

        //- Construct from stream (used by processor agglomeration gather)
        processorCyclicLduGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            Istream& is
        );


    //- Destructor
    virtual ~processorCyclicLduGAMGInterface() = default;
};

} // End namespace Foam

#endif

// ************************************************************************* //