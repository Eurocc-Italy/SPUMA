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
    Foam::turboPressureInletOutletVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This velocity inlet/outlet boundary condition is applied to velocity
    boundaries where the pressure is specified.  A zero-gradient condition is
    applied for outflow (as defined by the flux); for inflow, the velocity is
    obtained from the patch-face normal component of the internal-cell value.

    Direction components must be specified by inserting the required values
    in cylindrical coordinates.

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        axialComponent | Axial direction component | yes    |
        radialComponent | Radial direction component | yes  |
        tangentialComponent | tangential direction component | yes  |
        rotationAxis | Axis of rotation | no  |
        origin | origin of the cylindrical system of reference | no  |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            turboPressureInletOutletVelocity;
        phi             phi;
        // 45 deg in tangential-axial direction, no radial component
        // velocityComponents (radial tangential axial)
        axialComponent constant 0.7071;
        radialComponent constant 0.0;
        tangentialComponent constant 0.7071;
        rotationAxis (1 0 0);
        origin       (0 0 0);
        value        uniform (0 0 0);
    }
    \endverbatim

Note
    Sign conventions:
    - positive flux (out of domain): apply zero-gradient condition
    - negative flux (into of domain): derive from the flux in the patch-normal
      direction
    - sign of the rotation axis specifies the tangential direction
      according to right-hand-rule (thumb = radial, second finger = axial)

SourceFiles
    turboPressureInletOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef turboPressureInletOutletVelocityFvPatchVectorField_H
#define turboPressureInletOutletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class turboPressureInletOutletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class turboPressureInletOutletVelocityFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{
    // Private data

        //- Origin
        point origin_;

        //- Rotation axis;
        vector rotAxis_;

        //- Axial velocity
        autoPtr<Function1<scalar>> axialComponent_;

        //- Radial velocity
        autoPtr<Function1<scalar>> radialComponent_;

        //- Tangential velocity
        autoPtr<Function1<scalar>> tangentialComponent_;

        //- Flux field name
        word phiName_;

public:

    //- Runtime type information
    TypeName("turboPressureInletOutletVelocity");


    // Constructors

        //- Construct from patch and internal field
        turboPressureInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turboPressureInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turboPressureInletOutletVelocityFvPatchVectorField onto a new patch
        turboPressureInletOutletVelocityFvPatchVectorField
        (
            const turboPressureInletOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turboPressureInletOutletVelocityFvPatchVectorField
        (
            const turboPressureInletOutletVelocityFvPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        turboPressureInletOutletVelocityFvPatchVectorField
        (
            const turboPressureInletOutletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return fvPatchField<vector>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return fvPatchField<vector>::Clone(*this, iF);
        }


    // Member Functions

        //- True: this patch field is altered by assignment
        virtual bool assignable() const { return true; }


        // Access

            //- Return the name of phi
            const word& phiName() const
            {
                return phiName_;
            }

            //- Return reference to the name of phi to allow adjustment
            word& phiName()
            {
                return phiName_;
            }

            //- Reset the tangential velocity
            void updateVelocity();

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
