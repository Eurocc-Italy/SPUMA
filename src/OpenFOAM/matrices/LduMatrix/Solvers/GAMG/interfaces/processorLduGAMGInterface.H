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

Class
    Foam::processorLduGAMGInterface

Description
    Native GAMG agglomerated processor interface for templated LduMatrix GAMG.

SourceFiles
    processorLduGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_processorLduGAMGInterface_H
#define Foam_processorLduGAMGInterface_H

#include "LduGAMGInterface.H"
#include "processorLduInterface.H"

namespace Foam
{

class processorLduGAMGInterface
:
    public LduGAMGInterface,
    public processorLduInterface
{
    // Private data

        //- Communicator used for parallel communication.
        const label comm_;

        //- My processor rank in communicator.
        label myProcNo_;

        //- Neighbour processor rank in communicator.
        label neighbProcNo_;

        //- Face transformation tensor.
        tensorField forwardT_;

        //- Message tag.
        int tag_;


    // Private Member Functions

        processorLduGAMGInterface(const processorLduGAMGInterface&) = delete;
        void operator=(const processorLduGAMGInterface&) = delete;


public:

    TypeName("processorLduGAMGInterface");


    // Constructors

        processorLduGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& localRestrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );

        //- Construct from components (used by processor agglomeration combine)
        processorLduGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const labelUList& faceCells,
            const labelUList& faceRestrictAddressing,
            const label comm,
            const label myProcNo,
            const label neighbProcNo,
            const tensorField& forwardT,
            const int tag
        );

        //- Construct from stream (used by processor agglomeration gather)
        processorLduGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            Istream& is
        );


    //- Destructor
    virtual ~processorLduGAMGInterface() = default;


    // Member Functions

        // Interface transfer functions

            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;

            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData,
                const labelUList& faceCells
            ) const;

            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;


        // processorLduInterface implementation

            virtual label comm() const
            {
                return comm_;
            }

            virtual int myProcNo() const
            {
                return myProcNo_;
            }

            virtual int neighbProcNo() const
            {
                return neighbProcNo_;
            }

            virtual const tensorField& forwardT() const
            {
                return forwardT_;
            }

            virtual int tag() const
            {
                return tag_;
            }


        // I/O

            virtual void write(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //