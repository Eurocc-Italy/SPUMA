/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2022 Stefano Oliani
    Copyright (C) 2025 Cineca
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
    Foam::fvBlockMatrix

Description
    fvBlockMatrix is an extension of fvMatrix for block coupled types. It holds
    general insertion and retrieval tools for block coupling and specific
    functions for p-U coupling.

SourceFiles
    fvBlockMatrix.C
    fvBlockMatrixOperators.C
    fvBlockMatrixOperatorsTemplates.C
    fvBlockMatrixScalarScalarMatrix.C
    fvBlockMatrixScalarScalarMatrix.H
    fvBlockMatrixTensorVectorMatrix.C
    fvBlockMatrixTensorVectorMatrix.H
    fvBlockMatrixUpdateMatrixInterfaces.C


\*---------------------------------------------------------------------------*/

#ifndef Foam_fvBlockMatrix_H
#define Foam_fvBlockMatrix_H

#include "LduMatrix.H"
#include "tmp.H"
#include "fvMatricesFwd.H"
#include "surfaceFields.H"
#include "processorFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class sourceType, class blockType>
class fvBlockMatrix;

template<class Type2>
tmp<fvBlockMatrix<vector,Type2> > operator*
(
    const tmp<fvBlockMatrix<scalar,scalar> >& tA,
    const Type2& B
);


/*---------------------------------------------------------------------------*\
                           Class fvBlockMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class sourceType, class blockType>
class fvBlockMatrix
:
    public refCount,
    public LduMatrix<sourceType, blockType, blockType>
{
protected:

    // Protected Member Functions

        //- Add patch contribution to internal field
        template<class Type2>
        void addToInternalField
        (
            const labelUList& addr,
            const Field<Type2>& pf,
            Field<Type2>& intf
        ) const;

        template<class Type2>
        void addToInternalField
        (
            const labelUList& addr,
            const tmp<Field<Type2>>& tpf,
            Field<Type2>& intf
        ) const;

        //- Subtract patch contribution from internal field
        template<class Type2>
        void subtractFromInternalField
        (
            const labelUList& addr,
            const Field<Type2>& pf,
            Field<Type2>& intf
        ) const;

        template<class Type2>
        void subtractFromInternalField
        (
            const labelUList& addr,
            const tmp<Field<Type2>>& tpf,
            Field<Type2>& intf
        ) const;

        // Matrix completion functionality

        void addBoundaryDiag
        (
            fvMatrix<sourceType> matrix,
            scalarField& diag,
            const direction cmpt
        ) const;

        void addBoundarySource
        (
            fvMatrix<sourceType> matrix,
            Field<sourceType>& source,
            const bool couples=true
        ) const;

public:

    ClassName("fvBlockMatrix");

    // Constructors

        //- Construct given a field to solve for
        fvBlockMatrix(const lduMesh&);

        //- Construct as copy
        fvBlockMatrix(const fvBlockMatrix<sourceType,blockType>&);

        //- Clone
        tmp<fvBlockMatrix<sourceType,blockType>> clone() const;

    //- Destructor
    virtual ~fvBlockMatrix();

    // Member functions

        void insertBlock
        (
            const GeometricField<blockType, fvsPatchField, surfaceMesh>& leftField,
            const GeometricField<blockType, fvsPatchField, surfaceMesh>& rightField
        );

        void insertDissipationBlock
        (
            const GeometricField<blockType, fvsPatchField, surfaceMesh>& dissField
        );

        template<class psiType>
        void Amul
        (
            Field<sourceType>& Apsi,
            const GeometricField<psiType, fvPatchField, volMesh>& psi,
            const fvMesh& mesh
        ) const;

        template<class psiType>
        void AmulNoDiag
        (
            Field<sourceType>& Apsi,
            const GeometricField<psiType, fvPatchField, volMesh>& psi,
            const fvMesh& mesh
        ) const;

        //- Insert fvMatrix into this class (used if this is a
        // diagonal block without inter-component coupling
        // e.g. momentum equation)
        void insertEquation
        (
            fvMatrix<sourceType>& matrix
        );

        void insertMomentumEquation(fvVectorMatrix& matrix);

        void insertPressureEquation(fvScalarMatrix& matrix);

        void insertVelocityPressureCoupling
        (
            const fvBlockMatrix<vector, vector>& matrix
        );

        void insertPressureVelocityCoupling
        (
            const fvBlockMatrix<scalar, vector>& matrix
        );

        template<class psiType>
        void initProcessorInterfaces
        (
            const Field<psiType>& psiif,
            Field<psiType>& result,
            const LduInterfaceFieldPtrsList<psiType>& interfaces
        ) const;

        template<class psiType>
        void updateProcessorInterfaces
        (
            Field<sourceType>& result,
            const Pstream::commsTypes commsType,
            const processorFvPatchField<psiType>& procField,
            const fvMesh& mesh,
            const label patchi
        ) const;

    //- Operators

        void operator=(const fvBlockMatrix<sourceType, blockType>&);

        void operator=(const tmp<fvBlockMatrix<sourceType, blockType> >&);

        fvBlockMatrix<sourceType, blockType>& operator+=
        (
            const fvBlockMatrix<sourceType, blockType>&
        );

        fvBlockMatrix<sourceType, blockType>& operator+=
        (
            const tmp<fvBlockMatrix<sourceType, blockType> >&
        );

        fvBlockMatrix<sourceType, blockType>& operator-=
        (
            const fvBlockMatrix<sourceType, blockType>&
        );

        fvBlockMatrix<sourceType, blockType>& operator-=
        (
            const tmp<fvBlockMatrix<sourceType, blockType> >&
        );

    //- Friends

        template<class Type2>
        friend tmp<fvBlockMatrix<vector, Type2> > operator*
        (
            const tmp<fvBlockMatrix<scalar, scalar> >& tA, const Type2& B
        );
};


template<class sourceType, class blockType>
tmp<fvBlockMatrix<sourceType, blockType> > operator+
(
    const fvBlockMatrix<sourceType, blockType>&,
    const fvBlockMatrix<sourceType, blockType>&
);

template<class sourceType, class blockType>
tmp<fvBlockMatrix<sourceType, blockType> > operator+
(
    const tmp<fvBlockMatrix<sourceType, blockType> >&,
    const fvBlockMatrix<sourceType, blockType>&
);

template<class sourceType, class blockType>
tmp<fvBlockMatrix<sourceType, blockType> > operator+
(
    const fvBlockMatrix<sourceType, blockType>&,
    const tmp<fvBlockMatrix<sourceType, blockType> >&
);

template<class sourceType, class blockType>
tmp<fvBlockMatrix<sourceType, blockType> > operator+
(
    const tmp<fvBlockMatrix<sourceType, blockType> >&,
    const tmp<fvBlockMatrix<sourceType, blockType> >&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvBlockMatrix.C"
    #include "fvBlockMatrixUpdateMatrixInterfaces.C"
#endif

// Specialisation for scalars
#include "fvBlockMatrixScalarScalarMatrix.H"

// Specialisation for vectors
#include "fvBlockMatrixTensorVectorMatrix.H"

// Specialisation for vector4/tensor4 systems
#include "fvBlockMatrixVector4Tensor4Matrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
