/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "basicFvGeometryScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "primitiveMeshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicFvGeometryScheme, 0);
    addToRunTimeSelectionTable(fvGeometryScheme, basicFvGeometryScheme, dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicFvGeometryScheme::basicFvGeometryScheme
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvGeometryScheme(mesh, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::basicFvGeometryScheme::movePoints()
{
    fvGeometryScheme::movePoints();

    if (debug)
    {
        Pout<< "basicFvGeometryScheme::movePoints() : "
            << "recalculating primitiveMesh centres" << endl;
    }

    // Use lower level to calculate the geometry
    const_cast<fvMesh&>(mesh_).primitiveMesh::updateGeom();
}


Foam::tmp<Foam::surfaceScalarField> Foam::basicFvGeometryScheme::weights() const
{
    if (debug)
    {
        Pout<< "basicFvGeometryScheme::weights() : "
            << "Constructing weighting factors for face interpolation"
            << endl;
    }

    auto tweights =
        tmp<surfaceScalarField>::New
        (
            IOobject
            (
                "weights",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimless
        );

    auto& weights = tweights.ref();
    weights.setOriented();

    // Set local references to mesh data
    // Note that we should not use fvMesh sliced fields at this point yet
    // since this causes a loop when generating weighting factors in
    // coupledFvPatchField evaluation phase
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& C = mesh_.cellCentres();
    const vectorField& Sf = mesh_.faceAreas();

    // ... and reference to the internal field of the weighting factors
    scalarField& w = weights.primitiveFieldRef();

    foamExecutor exec;
    const auto SfPtr = Sf.cbegin();
    const auto CfPtr = Cf.cbegin();
    const auto CPtr = C.cbegin();
    const auto ownerPtr = owner.cbegin();
    const auto neighbourPtr = neighbour.cbegin();
    auto wPtr = w.begin();
    auto Lambda = [=](label facei)
    {
        // Note: mag in the dot-product.
        // For all valid meshes, the non-orthogonality will be less than
        // 90 deg and the dot-product will be positive.  For invalid
        // meshes (d & s <= 0), this will stabilise the calculation
        // but the result will be poor.
        scalar SfdOwn = mag(SfPtr[facei] & (CfPtr[facei] - CPtr[ownerPtr[facei]]));
        scalar SfdNei = mag(SfPtr[facei] & (CPtr[neighbourPtr[facei]] - CfPtr[facei]));

        if (mag(SfdOwn + SfdNei) > ROOTVSMALL)
        {
            wPtr[facei] = SfdNei/(SfdOwn + SfdNei);
        }
        else
        {
            wPtr[facei] = 0.5;
        }
    };
    exec.parallelFor(Lambda, owner.size());

    auto& wBf = weights.boundaryFieldRef();

    forAll(mesh_.boundary(), patchi)
    {
        mesh_.boundary()[patchi].makeWeights(wBf[patchi]);
    }

    if (debug)
    {
        Pout<< "basicFvGeometryScheme::weights : "
            << "Finished constructing weighting factors for face interpolation"
            << endl;
    }
    return tweights;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::basicFvGeometryScheme::deltaCoeffs() const
{
    if (debug)
    {
        Pout<< "basicFvGeometryScheme::deltaCoeffs() : "
            << "Constructing differencing factors array for face gradient"
            << endl;
    }

    // Force the construction of the weighting factors
    // needed to make sure deltaCoeffs are calculated for parallel runs.
    (void)mesh_.weights();

    auto tdeltaCoeffs =
        tmp<surfaceScalarField>::New
        (
            IOobject
            (
                "deltaCoeffs",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimless/dimLength
        );

    auto& deltaCoeffs = tdeltaCoeffs.ref();
    deltaCoeffs.setOriented();


    // Set local references to mesh data
    const vectorField& C = mesh_.cellCentres();
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    
    foamExecutor exec;
    const auto CPtr = C.cbegin();
    const auto ownerPtr = owner.cbegin();
    const auto neighbourPtr = neighbour.cbegin();
    auto deltaCoeffsPtr = deltaCoeffs.primitiveFieldRef().begin();
    auto Lambda = [=](label facei)
    {
        deltaCoeffsPtr[facei] = 1.0/mag(CPtr[neighbourPtr[facei]] - CPtr[ownerPtr[facei]]);
    };
    exec.parallelFor(Lambda, owner.size());

    auto& deltaCoeffsBf = deltaCoeffs.boundaryFieldRef();

    forAll(deltaCoeffsBf, patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];
        deltaCoeffsBf[patchi] = 1.0/mag(p.delta());

        // Optionally correct
        p.makeDeltaCoeffs(deltaCoeffsBf[patchi]);
    }

    return tdeltaCoeffs;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::basicFvGeometryScheme::nonOrthDeltaCoeffs() const
{
    if (debug)
    {
        Pout<< "basicFvGeometryScheme::nonOrthDeltaCoeffs() : "
            << "Constructing differencing factors array for face gradient"
            << endl;
    }

    // Force the construction of the weighting factors
    // needed to make sure deltaCoeffs are calculated for parallel runs.
    weights();

    auto tnonOrthDeltaCoeffs =
        tmp<surfaceScalarField>::New
        (
            IOobject
            (
                "nonOrthDeltaCoeffs",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimless/dimLength
        );

    auto& nonOrthDeltaCoeffs = tnonOrthDeltaCoeffs.ref();
    nonOrthDeltaCoeffs.setOriented();


    // Set local references to mesh data
    const volVectorField& C = mesh_.C();
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    foamExecutor exec;
    const auto CPtr = C.cbegin();
    const auto ownerPtr = owner.cbegin();
    const auto neighbourPtr = neighbour.cbegin();
    const auto SfPtr = Sf.cbegin();
    const auto magSfPtr = magSf.cbegin();
    auto nonOrthDeltaCoeffsPtr = nonOrthDeltaCoeffs.primitiveFieldRef().begin();
    auto Lambda = [=](label facei)
    {
        vector delta = CPtr[neighbourPtr[facei]] - CPtr[ownerPtr[facei]];
        vector unitArea = SfPtr[facei]/magSfPtr[facei];

        // Stabilised form for bad meshes
        nonOrthDeltaCoeffsPtr[facei] = 1.0/max(unitArea & delta, 0.05*mag(delta));
    };
    exec.parallelFor(Lambda, owner.size());

    auto& nonOrthDeltaCoeffsBf = nonOrthDeltaCoeffs.boundaryFieldRef();

    forAll(nonOrthDeltaCoeffsBf, patchi)
    {
        fvsPatchScalarField& patchDeltaCoeffs = nonOrthDeltaCoeffsBf[patchi];

        const fvPatch& p = patchDeltaCoeffs.patch();

        const vectorField patchDeltas(mesh_.boundary()[patchi].delta());

        const auto SfBfPtr = Sf.boundaryField()[patchi].cbegin();
        const auto magSfBfPtr = magSf.boundaryField()[patchi].cbegin();
        const auto patchDeltasPtr = patchDeltas.cbegin();
        auto patchDeltasCoeffsPtr = patchDeltaCoeffs.begin();
        auto Lambda = [=](label patchFacei)
        {
            vector unitArea = SfBfPtr[patchFacei]/magSfBfPtr[patchFacei];
            const vector& delta = patchDeltasPtr[patchFacei];

            patchDeltasCoeffsPtr[patchFacei] =
                1.0/max(unitArea & delta, 0.05*mag(delta));
        };
        exec.parallelFor(Lambda, p.size());


        // Optionally correct
        p.makeNonOrthoDeltaCoeffs(patchDeltaCoeffs);
    }
    return tnonOrthDeltaCoeffs;
}


Foam::tmp<Foam::surfaceVectorField>
Foam::basicFvGeometryScheme::nonOrthCorrectionVectors() const
{
    if (debug)
    {
        Pout<< "surfaceInterpolation::makeNonOrthCorrectionVectors() : "
            << "Constructing non-orthogonal correction vectors"
            << endl;
    }

    auto tnonOrthCorrectionVectors =
        tmp<surfaceVectorField>::New
        (
            IOobject
            (
                "nonOrthCorrectionVectors",
                mesh_.pointsInstance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimless
        );

    auto& corrVecs = tnonOrthCorrectionVectors.ref();
    corrVecs.setOriented();

    // Set local references to mesh data
    const volVectorField& C = mesh_.C();
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    tmp<surfaceScalarField> tNonOrthDeltaCoeffs(nonOrthDeltaCoeffs());
    const surfaceScalarField& NonOrthDeltaCoeffs = tNonOrthDeltaCoeffs();

    foamExecutor exec;
    const auto CPtr = C.cbegin();
    const auto ownerPtr = owner.cbegin();
    const auto neighbourPtr = neighbour.cbegin();
    const auto SfPtr = Sf.cbegin();
    const auto magSfPtr = magSf.cbegin();
    const auto NonOrthDeltaCoeffsPtr = NonOrthDeltaCoeffs.cbegin();
    auto corrVecsPtr = corrVecs.primitiveFieldRef().begin();
    auto Lambda = [=](label facei)
    {
        vector unitArea(SfPtr[facei]/magSfPtr[facei]);
        vector delta(CPtr[neighbourPtr[facei]] - CPtr[ownerPtr[facei]]);
        
        corrVecsPtr[facei] = unitArea - delta*NonOrthDeltaCoeffsPtr[facei];
    };
    exec.parallelFor(Lambda, owner.size());

    // Boundary correction vectors set to zero for boundary patches
    // and calculated consistently with internal corrections for
    // coupled patches

    auto& corrVecsBf = corrVecs.boundaryFieldRef();

    forAll(corrVecsBf, patchi)
    {
        fvsPatchVectorField& patchCorrVecs = corrVecsBf[patchi];

        const fvPatch& p = patchCorrVecs.patch();

        if (!patchCorrVecs.coupled())
        {
            patchCorrVecs = Zero;
        }
        else
        {
            const auto& patchNonOrthDeltaCoeffs =
                NonOrthDeltaCoeffs.boundaryField()[patchi];

            const vectorField patchDeltas(mesh_.boundary()[patchi].delta());

            const auto SfBfPtr = Sf.boundaryField()[patchi].cbegin();
            const auto magSfBfPtr = magSf.boundaryField()[patchi].cbegin();
            const auto patchDeltasPtr = patchDeltas.cbegin();
            const auto patchNonOrthDeltaCoeffsPtr = patchNonOrthDeltaCoeffs.cbegin();
            auto patchCorrVecsPtr = patchCorrVecs.begin();
            auto Lambda = [=](label patchFacei)
            {
                vector unitArea = SfBfPtr[patchFacei]/magSfBfPtr[patchFacei];
                const vector& delta = patchDeltasPtr[patchFacei];

                patchCorrVecsPtr[patchFacei] =
                    unitArea - delta*patchNonOrthDeltaCoeffsPtr[patchFacei];
            };
            exec.parallelFor(Lambda, p.size());
        }

        // Optionally correct
        p.makeNonOrthoCorrVectors(patchCorrVecs);
    }

    if (debug)
    {
        Pout<< "surfaceInterpolation::makeNonOrthCorrectionVectors() : "
            << "Finished constructing non-orthogonal correction vectors"
            << endl;
    }
    return tnonOrthCorrectionVectors;
}


bool Foam::basicFvGeometryScheme::updateGeom
(
    const pointField& points,
    const refPtr<pointField>& oldPoints,
    pointField& faceCentres,
    vectorField& faceAreas,
    pointField& cellCentres,
    scalarField& cellVolumes
) const
{
    primitiveMeshTools::makeFaceCentresAndAreas
    (
        mesh_,
        points,
        faceCentres,
        faceAreas
    );

    primitiveMeshTools::makeCellCentresAndVols
    (
        mesh_,
        faceCentres,
        faceAreas,
        cellCentres,
        cellVolumes
    );

    // Assume something has changed
    return true;
}


// ************************************************************************* //
