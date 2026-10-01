/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "turboPressureInletOutletVelocityFvPatchVectorField.H"
#include "diagTensor.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "vectorField.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turboPressureInletOutletVelocityFvPatchVectorField::
turboPressureInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    origin_(Zero),
    rotAxis_(Zero),
    axialComponent_(),
    radialComponent_(),
    tangentialComponent_(),
    phiName_("phi")
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::turboPressureInletOutletVelocityFvPatchVectorField::
turboPressureInletOutletVelocityFvPatchVectorField
(
    const turboPressureInletOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    rotAxis_(ptf.rotAxis_),
    axialComponent_(ptf.axialComponent_.clone()),
    radialComponent_(ptf.radialComponent_.clone()),
    tangentialComponent_(ptf.tangentialComponent_.clone()),
    phiName_(ptf.phiName_)
{}


Foam::turboPressureInletOutletVelocityFvPatchVectorField::
turboPressureInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    origin_(dict.lookup("origin")),
    rotAxis_(dict.lookup("rotationAxis")),
    axialComponent_(Function1<scalar>::New("axialComponent", dict, &db())),
    radialComponent_
    (
        Function1<scalar>::New("radialComponent", dict, &db())
    ),
    tangentialComponent_
    (
        Function1<scalar>::New("tangentialComponent", dict, &db())
    ),
    phiName_(dict.getOrDefault<word>("phi", "phi"))
{
    fvPatchFieldBase::readDict(dict);
    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::turboPressureInletOutletVelocityFvPatchVectorField::
turboPressureInletOutletVelocityFvPatchVectorField
(
    const turboPressureInletOutletVelocityFvPatchVectorField& pivpvf
)
:
    directionMixedFvPatchVectorField(pivpvf),
    origin_(pivpvf.origin_),
    rotAxis_(pivpvf.rotAxis_),
    axialComponent_(pivpvf.axialComponent_.clone()),
    radialComponent_(pivpvf.radialComponent_.clone()),
    tangentialComponent_(pivpvf.tangentialComponent_.clone()),
    phiName_(pivpvf.phiName_)
{}


Foam::turboPressureInletOutletVelocityFvPatchVectorField::
turboPressureInletOutletVelocityFvPatchVectorField
(
    const turboPressureInletOutletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(pivpvf, iF),
    origin_(pivpvf.origin_),
    rotAxis_(pivpvf.rotAxis_),
    axialComponent_(pivpvf.axialComponent_.clone()),
    radialComponent_(pivpvf.radialComponent_.clone()),
    tangentialComponent_(pivpvf.tangentialComponent_.clone()),
    phiName_(pivpvf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::turboPressureInletOutletVelocityFvPatchVectorField::
updateVelocity()
{
    const vectorField n(patch().nf());


    //Determining normal velocity from patchInternalField
    const scalarField normalVelocity =  this->patchInternalField() & rotAxis_;

    // Finding radial and tangential direction for each patch face center
    const symmTensor planeDir
    (
        1 - rotAxis_.x(), 
        0.0,            
        0.0,
        1 - rotAxis_.y(), 
        0.0,
        1 - rotAxis_.z()
    );
    vectorField radDir = transform(planeDir, patch().Cf()) - transform(planeDir, origin_);
    const scalarField magr(mag(radDir));
    radDir /= mag(radDir);

    const scalarField axialComp(axialComponent_->value(magr));
    const scalarField radComp(radialComponent_->value(magr));
    const scalarField tanComp(tangentialComponent_->value(magr));

    vectorField tanDir = rotAxis_ ^ radDir;

    // Calculating and assigning tangential velocity
    vectorField tangentialVelocity = (tanDir*tanComp + radDir*radComp)*normalVelocity/axialComp;
    refValue() = tangentialVelocity - n*(n & tangentialVelocity);
}

void Foam::turboPressureInletOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField>(phiName_);

    valueFraction() = neg(phip)*(I - sqr(patch().nf()));

    updateVelocity();
    directionMixedFvPatchVectorField::updateCoeffs();
    directionMixedFvPatchVectorField::evaluate();
}


void Foam::turboPressureInletOutletVelocityFvPatchVectorField::write
(
    Ostream& os
)
const
{
    fvPatchField<vector>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntry<vector>("rotationAxis", rotAxis_);
    os.writeEntry<point>("origin", origin_);
    axialComponent_->writeData(os);
    radialComponent_->writeData(os);
    tangentialComponent_->writeData(os);
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::turboPressureInletOutletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    tmp<vectorField> normalValue = transform(valueFraction(), refValue());
    tmp<vectorField> transformGradValue = transform(I - valueFraction(), pvf);
    fvPatchField<vector>::operator=(normalValue + transformGradValue);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        turboPressureInletOutletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
