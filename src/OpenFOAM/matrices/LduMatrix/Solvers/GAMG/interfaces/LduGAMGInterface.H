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
    Foam::LduGAMGInterface

Description
    Native coarse-interface base for templated LduMatrix GAMG.

SourceFiles
    LduGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_LduGAMGInterface_H
#define Foam_LduGAMGInterface_H

#include "autoPtr.H"
#include "lduInterfacePtrsList.H"
#include "primitiveFields.H"

namespace Foam
{

class Ostream;
class Istream;

class LduGAMGInterface
:
    public lduInterface
{
protected:

    // Protected data

        //- My index in the coarse interface list.
        const label index_;

        //- All interfaces for this coarse level.
        const lduInterfacePtrsList& coarseInterfaces_;

        //- Coarse-face to coarse-cell addressing.
        labelList faceCells_;

        //- Fine-face to coarse-face restriction addressing.
        labelList faceRestrictAddressing_;


    // Protected Member Functions

        //- No copy construct
        LduGAMGInterface(const LduGAMGInterface&) = delete;

        //- No copy assignment
        void operator=(const LduGAMGInterface&) = delete;


public:

    //- Runtime type information
    TypeName("LduGAMGInterface");


    // Selectors

        //- Construct a native coarse interface from a fine interface.
        static autoPtr<LduGAMGInterface> New
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& localRestrictAddressing,
            const labelField& neighbourRestrictAddressing,
            const label fineLevelIndex,
            const label coarseComm
        );

        //- Construct a native coarse interface from a stream (proc agglom).
        static autoPtr<LduGAMGInterface> New
        (
            const word& coupleType,
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            Istream& is
        );


    // Constructors

        //- Construct from interfaces, restrict addressing set later on.
        LduGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces
        )
        :
            index_(index),
            coarseInterfaces_(coarseInterfaces)
        {}


        //- Construct from interfaces and restrict addressing.
        LduGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            const labelUList& faceCells,
            const labelUList& faceRestrictAddressing
        )
        :
            index_(index),
            coarseInterfaces_(coarseInterfaces),
            faceCells_(faceCells),
            faceRestrictAddressing_(faceRestrictAddressing)
        {}


        //- Construct from interfaces and Istream (used by proc agglomeration).
        LduGAMGInterface
        (
            const label index,
            const lduInterfacePtrsList& coarseInterfaces,
            Istream& is
        )
        :
            index_(index),
            coarseInterfaces_(coarseInterfaces),
            faceCells_(is),
            faceRestrictAddressing_(is)
        {}


    //- Destructor
    virtual ~LduGAMGInterface() = default;


    // Member Functions

        // Access

            virtual label size() const
            {
                return faceCells_.size();
            }

            virtual label index() const
            {
                return index_;
            }

            virtual const lduInterfacePtrsList& coarseInterfaces() const
            {
                return coarseInterfaces_;
            }

            virtual const labelUList& faceCells() const
            {
                return faceCells_;
            }

            virtual const labelList& faceRestrictAddressing() const
            {
                return faceRestrictAddressing_;
            }

            virtual labelList& faceRestrictAddressing()
            {
                return faceRestrictAddressing_;
            }

            template<class Type>
            tmp<Field<Type>> interfaceInternalField
            (
                const UList<Type>& internalData
            ) const;

            template<class Type>
            tmp<Field<Type>> interfaceInternalField
            (
                const UList<Type>& internalData,
                const labelUList& faceCells
            ) const;

            template<class Type>
            void interfaceInternalField
            (
                const UList<Type>& internalData,
                List<Type>& result
            ) const;

            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData,
                const labelUList& faceCells
            ) const;


        // Agglomeration

            //- Merge the next level with this level.
            void combine(const LduGAMGInterface& coarseInterface);

            //- Agglomerate the fine coefficients for this interface.
            virtual tmp<scalarField> agglomerateCoeffs
            (
                const scalarField& fineCoeffs
            ) const;


        // I/O

            virtual void write(Ostream&) const = 0;
};

} // End namespace Foam

#ifdef NoRepository
    #include "LduGAMGInterfaceTemplates.C"
#endif

#endif

// ************************************************************************* //