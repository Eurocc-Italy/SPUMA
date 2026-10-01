/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turboTotalFvPatchScalarField.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turboTotalFvPatchScalarField::turboTotalFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    pName_("p"),
    p0_(p.size(), Zero),
    h0_(p.size(), Zero),
    sHPTable_("constant/sHPTable", bounds::normalBounding::CLAMP),
    pHSTable_("constant/pHSTable", bounds::normalBounding::CLAMP),
    THPTable_("constant/THPTable", bounds::normalBounding::CLAMP)
{}


Foam::turboTotalFvPatchScalarField::turboTotalFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, IOobjectOption::NO_READ),
    UName_(dict.getOrDefault<word>("U", "U")),
    pName_(dict.getOrDefault<word>("p", "p")),
    p0_("p0", dict, p.size()),
    h0_("h0", dict, p.size()),
    sHPTable_("constant/sHPTable", bounds::normalBounding::CLAMP),
    pHSTable_("constant/pHSTable", bounds::normalBounding::CLAMP),
    THPTable_("constant/THPTable", bounds::normalBounding::CLAMP)
{
    if (!this->readValueEntry(dict))
    {
        if (internalField().dimensions() == dimPressure)
        {
            fvPatchField<scalar>::operator=(p0_);
        }
        else if (internalField().dimensions() == dimTemperature)
        {
            scalarField T("T0",dict,p.size());
            fvPatchField<scalar>::operator=(T);
        }
        else
        {
            FatalErrorInFunction
                << " Incorrect field dimensions " << internalField().dimensions()
                << nl
                << "    Should be " << dimPressure
                << " for pressure boundary conditions" << nl
                << "    or " << dimTemperature
                << " for temperature boundary conditions," << nl
                << "    on patch " << this->patch().name()
                << " of field " << this->internalField().name()
                << " in file " << this->internalField().objectPath()
                << exit(FatalError);
        }
    }
}


Foam::turboTotalFvPatchScalarField::turboTotalFvPatchScalarField
(
    const turboTotalFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    pName_(ptf.pName_),
    p0_(ptf.p0_, mapper),
    h0_(ptf.h0_, mapper),
    sHPTable_("constant/sHPTable"),
    pHSTable_("constant/pHSTable"),
    THPTable_("constant/THPTable")
{}


Foam::turboTotalFvPatchScalarField::turboTotalFvPatchScalarField
(
    const turboTotalFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    pName_(tppsf.pName_),
    p0_(tppsf.p0_),
    h0_(tppsf.h0_),
    sHPTable_("constant/sHPTable"),
    pHSTable_("constant/pHSTable"),
    THPTable_("constant/THPTable")
{}


Foam::turboTotalFvPatchScalarField::turboTotalFvPatchScalarField
(
    const turboTotalFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    pName_(tppsf.pName_),
    p0_(tppsf.p0_),
    h0_(tppsf.h0_),
    sHPTable_("constant/sHPTable"),
    pHSTable_("constant/pHSTable"),
    THPTable_("constant/THPTable")
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turboTotalFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    p0_.autoMap(m);
    h0_.autoMap(m);
}


void Foam::turboTotalFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const turboTotalFvPatchScalarField& tiptf =
        refCast<const turboTotalFvPatchScalarField>(ptf);

    p0_.rmap(tiptf.p0_, addr);
    h0_.rmap(tiptf.h0_, addr);
}


void Foam::turboTotalFvPatchScalarField::updateCoeffs
(
    const scalarField& p0p,
    const scalarField& h0p,
    const vectorField& Up,
    const scalarField& pp
)
{
    if (updated())
    {
        return;
    }

    label patchSize = p0p.size();
    scalarField hp = h0p - 0.5*magSqr(Up);

    foamExecutor exec;
    auto& patchFld = *this;
    auto patchFldPtr = patchFld.begin();
    const auto p0pPtr = p0p.cbegin();
    const auto h0pPtr = h0p.cbegin();
    const auto UpPtr = Up.cbegin();
    const auto ppPtr = pp.cbegin();
    const auto hpPtr = hp.begin();

    if (internalField().dimensions() == dimPressure)
    {
        scalarField sp(patchSize);
        auto sHPTable(this->sHPTable_);
        auto pHSTable(this->pHSTable_);

        auto Lambda = [=](label facei)
        {
            scalar s0p = sHPTable(h0pPtr[facei], p0pPtr[facei]);
            // Assuming isentropic process inside the cell -> s0 = s
            patchFldPtr[facei] = pHSTable(hpPtr[facei], s0p);
        };
        exec.parallelFor(Lambda, patchFld.size());
    }
    else if (internalField().dimensions() == dimTemperature)
    {
        auto THPTable(this->THPTable_);

        auto Lambda = [=](label facei)
        {
            patchFldPtr[facei] = THPTable(hpPtr[facei], ppPtr[facei]);
        };
        exec.parallelFor(Lambda, patchFld.size());
    }
    else
    {
        FatalErrorInFunction
            << " Incorrect field dimensions " << internalField().dimensions()
            << nl
            << "    Should be " << dimPressure
            << " for pressure boundary conditions" << nl
            << "    or " << dimTemperature
            << " for temperature boundary conditions," << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::turboTotalFvPatchScalarField::updateCoeffs()
{
    updateCoeffs
    (
        p0(),
        h0(),
        patch().lookupPatchField<volVectorField>(UName()),
        patch().lookupPatchField<volScalarField>(pName())
    );
}


void Foam::turboTotalFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("p", "p", pName_);
    p0_.writeEntry("p0", os);
    h0_.writeEntry("h0", os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        turboTotalFvPatchScalarField
    );
}

// ************************************************************************* //
