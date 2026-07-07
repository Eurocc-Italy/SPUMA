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
    Foam::cyclicLduGAMGInterface

Description
    Native GAMG agglomerated cyclic interface for templated LduMatrix GAMG.

SourceFiles
    cyclicLduGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cyclicLduGAMGInterface_H
#define Foam_cyclicLduGAMGInterface_H

#include "LduGAMGInterface.H"
#include "cyclicLduInterface.H"

namespace Foam
{

class cyclicLduGAMGInterface
:
    public LduGAMGInterface,
    virtual public cyclicLduInterface
{
    // Private data

        //- Neighbour patch number.
        label neighbPatchID_;

        //- Ownership flag.
        bool owner_;

        //- Forward transformation tensor.
        tensorField forwardT_;

        //- Reverse transformation tensor.
        tensorField reverseT_;


    // Private Member Functions

        cyclicLduGAMGInterface(const cyclicLduGAMGInterface&) = delete;
        void operator=(const cyclicLduGAMGInterface&) = delete;


public:

    TypeName("cyclicLduGAMGInterface");


    // Constructors

        cyclicLduGAMGInterface
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
        cyclicLduGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const labelUList& faceCells,
            const labelUList& faceRestrictAddressing,
            const label neighbPatchID,
            const bool owner,
            const tensorField& forwardT,
            const tensorField& reverseT
        );

        //- Construct from stream (used by processor agglomeration gather)
        cyclicLduGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            Istream& is
        );


    //- Destructor
    virtual ~cyclicLduGAMGInterface() = default;


    // Member Functions

        // Interface transfer functions

            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;


        // cyclicLduInterface implementation

            virtual label neighbPatchID() const
            {
                return neighbPatchID_;
            }

            virtual bool owner() const
            {
                return owner_;
            }

            virtual const cyclicLduInterface& neighbPatch() const
            {
                return dynamic_cast<const cyclicLduGAMGInterface&>
                (
                    coarseInterfaces_[neighbPatchID_]
                );
            }

            virtual const tensorField& forwardT() const
            {
                return forwardT_;
            }

            virtual const tensorField& reverseT() const
            {
                return reverseT_;
            }


        // I/O

            virtual void write(Ostream& os) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //