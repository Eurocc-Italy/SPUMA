/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2025 CINECA (Ettore Fadiga)
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
    Foam::turboTotalFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition provides total boundary conditions for temperature
    or pressure fields.

Usage
    \table
        Property     | Description                | Required | Default value
        U            | Velocity field name        | no       | U
        p            | Pressure field name        | no       | phi
        h0           | Total enthalpy             | yes      |
        p0           | Total pressure             | yes      |
        T0           | Total temperature          | yes      |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            turboTotal;
        p0              uniform 1e5;
        h0              uniform 2e6;
        T0              uniform 300;
    }
    \endverbatim

    Requires s(H,p) table, p(H,s) table and T(H,p) table in the constant dir.
    Lookup tables must be named: sHPTable | pHSTable | THPTable

See also
    Foam::fixedValueFvPatchField

SourceFiles
    turboTotalFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turboTotalFvPatchScalarField_H
#define turboTotalFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class turboTotalFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turboTotalFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
        word UName_;

        //- Name of the pressure field
        word pName_;

        //- Total pressure
        scalarField p0_;

        //- Total enthalpy
        scalarField h0_;

        //- sHP table
        interpolation2DTable<scalar> sHPTable_;

        //- pHS table
        interpolation2DTable<scalar> pHSTable_;

        //- THP table
        interpolation2DTable<scalar> THPTable_;

public:

    //- Runtime type information
    TypeName("turboTotal");


    // Constructors

        //- Construct from patch and internal field
        turboTotalFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turboTotalFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turboTotalFvPatchScalarField
        //  onto a new patch
        turboTotalFvPatchScalarField
        (
            const turboTotalFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turboTotalFvPatchScalarField
        (
            const turboTotalFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        turboTotalFvPatchScalarField
        (
            const turboTotalFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member functions

        // Access

            //- Return the name of the velocity field
            const word& UName() const
            {
                return UName_;
            }

            //- Return reference to the name of the velocity field
            //  to allow adjustment
            word& UName()
            {
                return UName_;
            }

            //- Return the name of the pressure field
            const word& pName() const
            {
                return pName_;
            }

            //- Return reference to the name of the pressure field
            //  to allow adjustment
            word& pName()
            {
                return pName_;
            }

            //- Return the total pressure
            const scalarField& p0() const
            {
                return p0_;
            }

            //- Return reference to the total pressure to allow adjustment
            scalarField& p0()
            {
                return p0_;
            }

            //- Return the total enthalpy
            const scalarField& h0() const
            {
                return h0_;
            }

            //- Return reference to the total enthalpy to allow adjustment
            scalarField& h0()
            {
                return h0_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Inherit updateCoeffs from fixedValueFvPatchScalarField
            using fixedValueFvPatchScalarField::updateCoeffs;

            //- Update the coefficients associated with the patch field
            //  using the given patch total pressure and velocity fields
            virtual void updateCoeffs
            (
                const scalarField& p0p,
                const scalarField& h0p,
                const vectorField& Up,
                const scalarField& pp
            );

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
