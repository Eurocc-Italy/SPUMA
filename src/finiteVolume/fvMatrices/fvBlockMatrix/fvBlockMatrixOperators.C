/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2022 Stefano Oliani
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

#include "fvCFD.H"
#include "fvBlockMatrixOperators.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvm
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<fvBlockMatrix<vector,vector>> grad
(
    const volScalarField& vf,
    const word& name
)
{
    tmp<surfaceInterpolationScheme<scalar>> tInterpolationScheme
    (
        surfaceInterpolationScheme<scalar>::New
        (
            vf.mesh(),
            vf.mesh().interpolationScheme(name)
        )
    );

    tmp<surfaceScalarField> tweights = tInterpolationScheme().weights(vf);
    const scalarField& wIn = tweights().internalField();

    const fvMesh& mesh = vf.mesh();

    typedef fvBlockMatrix<vector,vector> blockCoupledMatrix;

    tmp<blockCoupledMatrix> tbs
    (
        new blockCoupledMatrix(mesh)
    );

    fvBlockMatrix<vector,vector>& bs = tbs.ref();

    vectorField& source = bs.source();

    // Grab ldu parts of block matrix
    vectorField& d = bs.diag();
    vectorField& u = bs.upper();
    vectorField& l = bs.lower();

    const vectorField& SfIn = mesh.Sf().internalField();

    l = -wIn*SfIn;
    u = l + SfIn;
    bs.negSumDiag();

    bs.interfacesUpper().resize(mesh.boundary().size());
    bs.interfacesLower().resize(mesh.boundary().size());

    foamExecutor exec;

    // Boundary contributions
    forAll (vf.boundaryField(), patchI)
    {
        const fvPatchScalarField& pf = vf.boundaryField()[patchI];
        const fvPatch& patch = pf.patch();
        const vectorField& pSf = patch.Sf();
        const fvsPatchScalarField& pw = tweights().boundaryField()[patchI];
        const labelList& fc = patch.faceCells();

        const vectorField pcl = -pw*pSf;
        const vectorField pcu = pcl + pSf;

        bs.interfacesUpper().set(patchI, new vectorField(pcu));
        bs.interfacesLower().set(patchI, new vectorField(pcl));

        const scalarField internalCoeffs(pf.valueInternalCoeffs(pw));

        vector* __restrict__ dPtr = d.begin();
        const label* const __restrict__ fcPtr = fc.cbegin();
        const vector* const __restrict__ pSfPtr = pSf.cbegin();
        const scalar* const __restrict__ internalCoeffsPtr = internalCoeffs.cbegin();

        // Diag contribution
        auto LambdaDiag = [=](label faceI)
        {
            foamAtomic::AtomicAdd(dPtr[fcPtr[faceI]], internalCoeffsPtr[faceI]*pSfPtr[faceI]);
        };
        exec.parallelFor(LambdaDiag, pf.size());

        if (patch.coupled())
        {
        }
        else
        {
            const scalarField boundaryCoeffs(pf.valueBoundaryCoeffs(pw));

            const scalar* const __restrict__ boundaryCoeffsPtr = boundaryCoeffs.cbegin();
            vector* __restrict__ sourcePtr = source.begin();

            // Boundary contribution
            auto LambdaSrc = [=](label faceI)
            {
                foamAtomic::AtomicAdd(sourcePtr[fcPtr[faceI]], -boundaryCoeffsPtr[faceI]*pSfPtr[faceI]);
            };
            exec.parallelFor(LambdaSrc, pf.size());
        }
    }

    // Interpolation schemes with corrections not accounted for

    return tbs;
}

tmp<fvBlockMatrix<vector,vector>> grad
(
    const volScalarField& vf
)
{
    return fvm::grad
    (
        vf,
        "grad(" + vf.name() + ')'
    );
}

tmp<fvBlockMatrix<scalar,vector>> div
(
    const volVectorField& vf
)
{
    const fvMesh& mesh = vf.mesh();
    const surfaceScalarField& meshWeights = mesh.weights();

    tmp<fvBlockMatrix<scalar, vector>> tbs
    (
        new fvBlockMatrix<scalar, vector>(mesh)
    );

    fvBlockMatrix<scalar, vector>& bs = tbs.ref();

    scalarField& source = bs.source();

    // Grab ldu parts of block matrix
    vectorField& d = bs.diag();
    vectorField& u = bs.upper();
    vectorField& l = bs.lower();

    const vectorField& SfIn = mesh.Sf().internalField();
    const scalarField& wIn = meshWeights.internalField();

    l = -wIn*SfIn;
    u = l + SfIn;
    bs.negSumDiag();

    bs.interfacesUpper().resize(mesh.boundary().size());
    bs.interfacesLower().resize(mesh.boundary().size());

    foamExecutor exec;

    // Boundary contributions
    forAll (vf.boundaryField(), patchI)
    {
        const fvPatchVectorField& pf = vf.boundaryField()[patchI];
        const fvPatch& patch = pf.patch();
        const vectorField& pSf = patch.Sf();
        const fvsPatchScalarField& pw = meshWeights.boundaryField()[patchI];
        const labelList& fc = patch.faceCells();

        const vectorField pcl = -pw*pSf;
        const vectorField pcu = pcl + pSf;

        bs.interfacesUpper().set(patchI, new vectorField(pcu));
        bs.interfacesLower().set(patchI, new vectorField(pcl));

        const vectorField internalCoeffs(pf.valueInternalCoeffs(pw));

        vector* __restrict__ dPtr = d.begin();
        const label* const __restrict__ fcPtr = fc.cbegin();
        const vector* const __restrict__ pSfPtr = pSf.cbegin();
        const vector* const __restrict__ internalCoeffsPtr = internalCoeffs.cbegin();

        // Diag contribution
        auto LambdaDiag = [=](label faceI)
        {
            foamAtomic::AtomicAdd(dPtr[fcPtr[faceI]], cmptMultiply(internalCoeffsPtr[faceI], pSfPtr[faceI]));
        };
        exec.parallelFor(LambdaDiag, pf.size());

        if (patch.coupled())
        {
        }
        else
        {
            const vectorField boundaryCoeffs(pf.valueBoundaryCoeffs(pw));

            const vector* const __restrict__ boundaryCoeffsPtr = boundaryCoeffs.cbegin();
            scalar* __restrict__ sourcePtr = source.begin();

            // Boundary contribution
            auto LambdaSrc = [=](label faceI)
            {
                foamAtomic::AtomicAdd(sourcePtr[fcPtr[faceI]], -boundaryCoeffsPtr[faceI] & pSfPtr[faceI]);
            };
            exec.parallelFor(LambdaSrc, pf.size());
        }
    }

    // Interpolation schemes with corrections not accounted for

    return tbs;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace fvj
{

tmp<fvBlockMatrix<scalar, vector> >
grad
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const geometricOneField&
)
{
    const fvMesh& mesh = w.mesh();

    tmp<fvBlockMatrix<scalar, vector > > tmx
    (
        new fvBlockMatrix<scalar, vector>(mesh)
    );
    fvBlockMatrix<scalar, vector>& mx = tmx.ref();

    tmp< surfaceVectorField > lsf = -w*mesh.Sf();
    tmp< surfaceVectorField > usf = mesh.Sf() + lsf();

    Field<vector>& upp = mx.upper();
    Field<vector>& low = mx.lower();
    Field<vector>& diag = mx.diag();

    forAll(upp, facei)
    {
        upp[facei] = usf()[facei];
        low[facei] = lsf()[facei];
    }

    mx.negSumDiag();

    mx.interfacesUpper().resize(mesh.boundary().size());
    mx.interfacesLower().resize(mesh.boundary().size());

    forAll(mesh.boundary(), patchi)
    {
        mx.interfacesUpper().set(patchi, new vectorField(usf->boundaryField()[patchi]));
        mx.interfacesLower().set(patchi, new vectorField(lsf->boundaryField()[patchi]));

        // Don't include physical boundaries because those are dealt with separately
        if (mesh.boundary()[patchi].coupled())
        {

            Field<vector>& low =
                mx.interfacesLower()[patchi];

            const labelUList& bOwner =
                mesh.boundary()[patchi].faceCells();

            forAll(bOwner, facei)
            {
                const label own = bOwner[facei];
                diag[own] -= low[facei];
            }
        }
    }
    return tmx;
}

} // End namespace fvj

} // End namespace Foam

// ************************************************************************* //
