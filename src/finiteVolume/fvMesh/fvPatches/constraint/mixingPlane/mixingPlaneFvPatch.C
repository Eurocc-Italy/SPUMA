/*---------------------------------------------------------------------------*\

License
    This file is part of the TurboWG project:
    https://sourceforge.net/projects/turbowg

    The files from the TurboWG project are free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    This source code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this source code.  If not, see <http://www.gnu.org/licenses/>.

    Special notice:
    This file is part of a port of the mixingPlane interface source code from
    the foam-extend project (https://sourceforge.net/projects/foam-extend/) to
    OpenFOAM (https://www.openfoam.com/). Most of the original source code has
    been preserved, including the comments and authorship/IP notices.

Contributors/Copyright:
    Martin Beaudoin, (2023)

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    MixingPlane patch, providing coupling between arbitrary patches which
    belong to the same fvMesh using a mixing plane averaging technique

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "mixingPlaneFvPatch.H"
#include "fvPatchFields.H"
#include "fvsPatchFields.H"
#include "fvBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mixingPlaneFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, mixingPlaneFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Make patch weighting factors
void Foam::mixingPlaneFvPatch::makeWeights(scalarField& w) const
{
    // Calculation of weighting factors is performed from the master
    // position, using reconstructed shadow cell centres
    if (mixingPlanePolyPatch_.master())
    {
        const vectorField n(nf());

        // Note: mag in the dot-product.
        // For all valid meshes, the non-orthogonality will be less than
        // 90 deg and the dot-product will be positive.  For invalid
        // meshes (d & s <= 0), this will stabilise the calculation
        // but the result will be poor.  HJ, 24/Aug/2011
        const scalarField nfc
        (
            mag(n & (mixingPlanePolyPatch_.reconFaceCellCentres() - Cf()))
        );

        w = nfc/(mag(n & (Cf() - Cn())) + nfc);
    }
    else
    {
        // Pick up weights from the master side
        scalarField masterWeights(shadow().size());

        shadow().makeWeights(masterWeights);

        const scalarField oneMinusW(1 - masterWeights);

        w = interpolate(oneMinusW);
    }
}


void Foam::mixingPlaneFvPatch::makeDeltaCoeffs(scalarField& dc) const
{
    if (mixingPlanePolyPatch_.master())
    {
        // Stabilised form for bad meshes.  HJ, 24/Aug/2011
        const vectorField d(delta());

        dc = 1.0/max(nf() & d, 0.05*mag(d));
    }
    else
    {
        scalarField masterDeltas(shadow().size());

        shadow().makeDeltaCoeffs(masterDeltas);

        dc = interpolate(masterDeltas);
    }
}


void Foam::mixingPlaneFvPatch::makeCorrVecs(vectorField& cv) const
{
    cv = vector::zero;
#if 0
    // Full non-orthogonality treatment

    // Calculate correction vectors on coupled patches
    const scalarField& patchDeltaCoeffs = deltaCoeffs();

    const vectorField patchDeltas(delta());
    const vectorField n(nf());
    cv = n - patchDeltas*patchDeltaCoeffs;
#endif
}


const Foam::mixingPlaneFvPatch& Foam::mixingPlaneFvPatch::shadow() const
{
    const fvPatch& p =
        this->boundaryMesh()[mixingPlanePolyPatch_.shadowIndex()];

    return refCast<const mixingPlaneFvPatch>(p);
}


// Return delta (P to N) vectors across coupled patch
Foam::tmp<Foam::vectorField> Foam::mixingPlaneFvPatch::delta() const
{
    if (mixingPlanePolyPatch_.master())
    {
        return mixingPlanePolyPatch_.reconFaceCellCentres() - Cn();
    }
    else
    {
        tmp<vectorField> tDelta = interpolate
        (
            shadow().Cn()
          - mixingPlanePolyPatch_.shadow().reconFaceCellCentres()
        );

        return tDelta;
    }
}


bool Foam::mixingPlaneFvPatch::master() const
{
    return mixingPlanePolyPatch_.master();
}


Foam::label Foam::mixingPlaneFvPatch::shadowIndex() const
{
    return mixingPlanePolyPatch_.shadowIndex();
}


const Foam::mixingPlaneLduInterface& Foam::mixingPlaneFvPatch::shadowInterface() const
{
    const fvPatch& p =
        this->boundaryMesh()[mixingPlanePolyPatch_.shadowIndex()];

    return refCast<const mixingPlaneLduInterface>(p);
}


const Foam::labelListList& Foam::mixingPlaneFvPatch::ggiAddressing() const
{
    if (mixingPlanePolyPatch_.master())
    {
        return mixingPlanePolyPatch_.patchToPatch().masterPatchToProfileAddr();
    }
    else
    {
        return mixingPlanePolyPatch_.patchToPatch().slavePatchToProfileAddr();
    }
}


const Foam::scalarListList& Foam::mixingPlaneFvPatch::ggiWeights() const
{
    if (mixingPlanePolyPatch_.master())
    {
        return mixingPlanePolyPatch_.patchToPatch().
            masterPatchToProfileWeights();
    }
    else
    {
        return mixingPlanePolyPatch_.patchToPatch().
            slavePatchToProfileWeights();
    }
}


Foam::tmp<Foam::labelField> Foam::mixingPlaneFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}

Foam::tmp<Foam::labelField> Foam::mixingPlaneFvPatch::interfaceInternalField
(
    const labelUList& internalData,
    const labelUList& faceCells
) const
{
#if OPENFOAM < 2302
    return patchInternalField(internalData, faceCells);
#else
    auto tpfld = tmp<labelField>::New();
    patchInternalField(internalData, faceCells, tpfld.ref());
    return tpfld;
#endif
}

void Foam::mixingPlaneFvPatch::initTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& interfaceData
) const
{
    labelTransferBuffer_ = interfaceData;
}


Foam::tmp<Foam::labelField> Foam::mixingPlaneFvPatch::transfer
(
    const Pstream::commsTypes,
    const labelUList& interfaceData
) const
{
    return this->shadow().labelTransferBuffer();
}


void Foam::mixingPlaneFvPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    labelUList& iF
) const
{
    labelTransferBuffer_ = patchInternalField(iF);
}


Foam::tmp<Foam::labelField> Foam::mixingPlaneFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes,
    const labelUList& iF
) const
{
    return shadow().labelTransferBuffer();
}


// ************************************************************************* //
