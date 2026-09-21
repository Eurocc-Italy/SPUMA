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

#ifndef Foam_fvBlockMatrixOperators_H
#define Foam_fvBlockMatrixOperators_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "surfaceInterpolationScheme.H"
#include "fvBlockMatrix.H"
#include "tmp.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvm functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvm
{
    tmp<fvBlockMatrix<vector,vector>> grad
    (
        const volScalarField& vf,
        const word& name
    );

    tmp<fvBlockMatrix<vector,vector>> grad
    (
        const volScalarField& vf
    );

    tmp<fvBlockMatrix<scalar,vector>> div
    (
        const volVectorField& vf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Namespace of operators for density-based solvers
namespace fvj
{

//- Divergence of a cell Jacobian field interpolated to faces with
//  the specified weighting

template<class sourceType, class blockType>
tmp
<
    fvBlockMatrix
    <
        sourceType,
        typename innerProduct<blockType,vector>::type
    >
>
div
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const GeometricField<blockType, fvPatchField, volMesh>& vf,
    sourceType dummyValue  //Needed only to decide the source type
);

template<class sourceType, class blockType>
tmp
<
    fvBlockMatrix
    <
        sourceType,
        typename innerProduct<blockType,vector>::type
    >
>
div
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const tmp<GeometricField<blockType, fvPatchField, volMesh> >& tvf,
    sourceType dummyValue  //Needed only to decide the source type
);

template<class blockType>
tmp<fvBlockMatrix<blockType, blockType> >
div
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const GeometricField<blockType, fvsPatchField, surfaceMesh>& sf
);

template<class blockType>
tmp<fvBlockMatrix<blockType, blockType> >
div
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const tmp<GeometricField<blockType, fvsPatchField, surfaceMesh> >& tsf
);

template<class sourceType, class blockType>
tmp<fvBlockMatrix<sourceType, typename outerProduct<blockType,vector>::type> >
grad
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const GeometricField<blockType, fvPatchField, volMesh>& vf,
    sourceType dummyValue  //Needed only to decide the source type
);

template<class sourceType, class blockType>
tmp<fvBlockMatrix<sourceType, typename outerProduct<blockType,vector>::type> >
grad
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const tmp<GeometricField<blockType, fvPatchField, volMesh> >& tvf,
    sourceType dummyValue  //Needed only to decide the source type
);

tmp<fvBlockMatrix<scalar, vector> >
grad
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const geometricOneField&
);

template<class blockType>
tmp<fvBlockMatrix<blockType, typename outerProduct<blockType,vector>::type> >
gradT
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const GeometricField<blockType, fvPatchField, volMesh>& vf
);

template<class blockType>
tmp<fvBlockMatrix<blockType, typename outerProduct<blockType,vector>::type> >
gradT
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const tmp<GeometricField<blockType, fvPatchField, volMesh> >& tvf
);

template<class blockType>
tmp<fvBlockMatrix<blockType, blockType> >
laplacian
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& sf,
    const GeometricField<blockType, fvPatchField, volMesh>& vf,
    bool includePhysicalBoundaries=false
);

template<class blockType>
tmp<fvBlockMatrix<blockType, blockType> >
laplacian
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& sf,
    const tmp<GeometricField<blockType, fvPatchField, volMesh> >& tvf,
    bool includePhysicalBoundaries=false
);

template<class blockType>
tmp<fvBlockMatrix<blockType, blockType> >
laplacian
(
    const tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> >& tsf,
    const GeometricField<blockType, fvPatchField, volMesh>& vf,
    bool includePhysicalBoundaries=false
);

template<class blockType>
tmp<fvBlockMatrix<blockType, blockType> >
laplacian
(
    const tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> >& tsf,
    const tmp<GeometricField<blockType, fvPatchField, volMesh> >& tvf,
    bool includePhysicalBoundaries=false
);

template<class blockType>
tmp<fvBlockMatrix<blockType, blockType> >
laplacian
(
    const GeometricField<blockType, fvsPatchField, surfaceMesh>& sf,
    const geometricOneField&,
    bool includePhysicalBoundaries=false
);

template<class blockType>
tmp<fvBlockMatrix<blockType, blockType> >
laplacian
(
    const tmp<GeometricField<blockType, fvsPatchField, surfaceMesh> >& tsf,
    const geometricOneField&,
    bool includePhysicalBoundaries=false
);

template<class sourceType, class blockType>
tmp<fvBlockMatrix<sourceType, blockType> >
laplacian
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& sf,
    const GeometricField<blockType, fvPatchField, volMesh>& tvf,
    sourceType dummyValue  //Needed only to decide the source type
);

template<class sourceType, class blockType>
tmp<fvBlockMatrix<sourceType, blockType> >
laplacian
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& sf,
    const tmp<GeometricField<blockType, fvPatchField, volMesh> >& tvf,
    sourceType dummyValue  //Needed only to decide the source type
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvj

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvBlockMatrixOperatorsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
