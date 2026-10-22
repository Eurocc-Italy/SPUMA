/*---------------------------------------------------------------------------*\
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
    Foam::Tensor4

Description
    A templated (4 x 4) tensor derived from MatrixSpace.

    Supports component-wise VectorSpace operations, transpose, Tensor4 inner
    product and Tensor4/Vector4 inner product. Determinant, inverse and
    symmetric/spherical forms are not implemented.

SourceFiles
    Tensor4I.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_Tensor4_H
#define Foam_Tensor4_H

#include "deviceM.H"
#include "MatrixSpace.H"
#include "Vector4.H"
#include "tensor.H"

namespace Foam
{

template<class Cmpt>
class Tensor4
:
    public MatrixSpace<Tensor4<Cmpt>, Cmpt, 4, 4>
{
public:

    // Typedefs

        //- Equivalent type of labels used for valid component indexing
    typedef Tensor4<label> labelType;


    // Member Constants

        //- Rank of Tensor4 is 2
    static constexpr direction rank = 2;


    //- Component labeling enumeration
    enum components
    {
        XX, XY, XZ, XK,
        YX, YY, YZ, YK,
        ZX, ZY, ZZ, ZK,
        KX, KY, KZ, KK
    };


    // Generated Methods

        //- Default construct
    FOAM_DEVICE Tensor4() = default;

        //- Copy construct
    FOAM_DEVICE Tensor4(const Tensor4&) = default;

        //- Copy assignment
    FOAM_DEVICE Tensor4& operator=(const Tensor4&) = default;


    // Constructors

        //- Construct initialized to zero
    FOAM_DEVICE inline Tensor4(Foam::zero);

        //- Construct from a MatrixSpace with the same dimensions
    template<class Cmpt2>
    FOAM_DEVICE inline Tensor4
    (
        const MatrixSpace<Tensor4<Cmpt2>, Cmpt2, 4, 4>& matrix
    );

        //- Construct from an equivalent VectorSpace
    template<class Cmpt2>
    FOAM_DEVICE inline Tensor4
    (
        const VectorSpace<Tensor4<Cmpt2>, Cmpt2, 16>& vectorSpace
    );

        //- Construct from the sixteen components
    FOAM_DEVICE inline Tensor4
    (
        const Cmpt& txx, const Cmpt& txy, const Cmpt& txz, const Cmpt& txk,
        const Cmpt& tyx, const Cmpt& tyy, const Cmpt& tyz, const Cmpt& tyk,
        const Cmpt& tzx, const Cmpt& tzy, const Cmpt& tzz, const Cmpt& tzk,
        const Cmpt& tkx, const Cmpt& tky, const Cmpt& tkz, const Cmpt& tkk
    );

        //- Construct from Istream
    inline explicit Tensor4(Istream& is);


    // Member Functions

    // Component Access

        //- Access the xx component
    FOAM_DEVICE const Cmpt& xx() const noexcept { return this->v_[XX]; }

        //- Access the xy component
    FOAM_DEVICE const Cmpt& xy() const noexcept { return this->v_[XY]; }

        //- Access the xz component
    FOAM_DEVICE const Cmpt& xz() const noexcept { return this->v_[XZ]; }

        //- Access the xk component
    FOAM_DEVICE const Cmpt& xk() const noexcept { return this->v_[XK]; }

        //- Access the yx component
    FOAM_DEVICE const Cmpt& yx() const noexcept { return this->v_[YX]; }

        //- Access the yy component
    FOAM_DEVICE const Cmpt& yy() const noexcept { return this->v_[YY]; }

        //- Access the yz component
    FOAM_DEVICE const Cmpt& yz() const noexcept { return this->v_[YZ]; }

        //- Access the yk component
    FOAM_DEVICE const Cmpt& yk() const noexcept { return this->v_[YK]; }

        //- Access the zx component
    FOAM_DEVICE const Cmpt& zx() const noexcept { return this->v_[ZX]; }

        //- Access the zy component
    FOAM_DEVICE const Cmpt& zy() const noexcept { return this->v_[ZY]; }

        //- Access the zz component
    FOAM_DEVICE const Cmpt& zz() const noexcept { return this->v_[ZZ]; }

        //- Access the zk component
    FOAM_DEVICE const Cmpt& zk() const noexcept { return this->v_[ZK]; }

        //- Access the kx component
    FOAM_DEVICE const Cmpt& kx() const noexcept { return this->v_[KX]; }

        //- Access the ky component
    FOAM_DEVICE const Cmpt& ky() const noexcept { return this->v_[KY]; }

        //- Access the kz component
    FOAM_DEVICE const Cmpt& kz() const noexcept { return this->v_[KZ]; }

        //- Access the kk component
    FOAM_DEVICE const Cmpt& kk() const noexcept { return this->v_[KK]; }

        //- Access the xx component
    FOAM_DEVICE Cmpt& xx() noexcept { return this->v_[XX]; }

        //- Access the xy component
    FOAM_DEVICE Cmpt& xy() noexcept { return this->v_[XY]; }

        //- Access the xz component
    FOAM_DEVICE Cmpt& xz() noexcept { return this->v_[XZ]; }

        //- Access the xk component
    FOAM_DEVICE Cmpt& xk() noexcept { return this->v_[XK]; }

        //- Access the yx component
    FOAM_DEVICE Cmpt& yx() noexcept { return this->v_[YX]; }

        //- Access the yy component
    FOAM_DEVICE Cmpt& yy() noexcept { return this->v_[YY]; }

        //- Access the yz component
    FOAM_DEVICE Cmpt& yz() noexcept { return this->v_[YZ]; }

        //- Access the yk component
    FOAM_DEVICE Cmpt& yk() noexcept { return this->v_[YK]; }

        //- Access the zx component
    FOAM_DEVICE Cmpt& zx() noexcept { return this->v_[ZX]; }

        //- Access the zy component
    FOAM_DEVICE Cmpt& zy() noexcept { return this->v_[ZY]; }

        //- Access the zz component
    FOAM_DEVICE Cmpt& zz() noexcept { return this->v_[ZZ]; }

        //- Access the zk component
    FOAM_DEVICE Cmpt& zk() noexcept { return this->v_[ZK]; }

        //- Access the kx component
    FOAM_DEVICE Cmpt& kx() noexcept { return this->v_[KX]; }

        //- Access the ky component
    FOAM_DEVICE Cmpt& ky() noexcept { return this->v_[KY]; }

        //- Access the kz component
    FOAM_DEVICE Cmpt& kz() noexcept { return this->v_[KZ]; }

        //- Access the kk component
    FOAM_DEVICE Cmpt& kk() noexcept { return this->v_[KK]; }


    // Member Operators

        //- Inherit MatrixSpace assignment operators
    using Tensor4::msType::operator=;

        //- Assign from an equivalent vector space
    template<class Cmpt2>
    FOAM_DEVICE inline void operator=
    (
        const VectorSpace<Tensor4<Cmpt2>, Cmpt2, 16>& vectorSpace
    );
};


//- Return the inner product of a Tensor4 and a Vector4
template<class Cmpt>
FOAM_DEVICE inline Vector4<Cmpt> dot
(
    const Tensor4<Cmpt>& matrix,
    const Vector4<Cmpt>& vector
)
{
    return matrix & vector;
}


//- Transform the velocity components of a Vector4; preserve the scalar K
template<class Cmpt>
FOAM_DEVICE inline Vector4<Cmpt> transform
(
    const tensor& tt,
    const Vector4<Cmpt>& vector
)
{
    return Vector4<Cmpt>
    (
        tt.xx()*vector.x() + tt.xy()*vector.y() + tt.xz()*vector.z(),
        tt.yx()*vector.x() + tt.yy()*vector.y() + tt.yz()*vector.z(),
        tt.zx()*vector.x() + tt.zy()*vector.y() + tt.zz()*vector.z(),
        vector.k()
    );
}


//- Solve a Tensor4 block system with partial pivoting
template<class Cmpt>
FOAM_DEVICE inline Vector4<Cmpt> solve
(
    const Tensor4<Cmpt>& matrix,
    const Vector4<Cmpt>& rhs
);


//- Return the inverse of a Tensor4, column by column via solve()
template<class Cmpt>
FOAM_DEVICE inline Tensor4<Cmpt> inv(const Tensor4<Cmpt>& matrix);


// * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<Tensor4<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous labels if component type is a contiguous label
template<class Cmpt>
struct is_contiguous_label<Tensor4<Cmpt>> : is_contiguous_label<Cmpt> {};

//- Data are contiguous scalars if component type is a contiguous scalar
template<class Cmpt>
struct is_contiguous_scalar<Tensor4<Cmpt>> : is_contiguous_scalar<Cmpt> {};

//- The transpose of a Tensor4 is a Tensor4
template<class Cmpt>
class typeOfTranspose<Cmpt, Tensor4<Cmpt>>
{
public:
    typedef Tensor4<Cmpt> type;
};

//- The inner product of two Tensor4 objects is a Tensor4
template<class Cmpt>
class typeOfInnerProduct<Cmpt, Tensor4<Cmpt>, Tensor4<Cmpt>>
{
public:
    typedef Tensor4<Cmpt> type;
};

//- The inner product of Tensor4 and Vector4 is a Vector4
template<class Cmpt>
class typeOfInnerProduct<Cmpt, Tensor4<Cmpt>, Vector4<Cmpt>>
{
public:
    typedef Vector4<Cmpt> type;
};

//- Use solveScalar components when solving for a Tensor4
template<class Cmpt>
class typeOfSolve<Tensor4<Cmpt>>
{
public:
    typedef Tensor4<solveScalar> type;
};


//- Transform the spatial dimensions, leaving the fourth axis invariant
template<class Cmpt>
FOAM_DEVICE_I Tensor4<Cmpt> transform
(
    const tensor& rotation,
    const Tensor4<Cmpt>& value
)
{
    Tensor4<Cmpt> result(Zero);

    for (direction row = 0; row < 4; ++row)
    {
        for (direction column = 0; column < 4; ++column)
        {
            for (direction innerRow = 0; innerRow < 4; ++innerRow)
            {
                const scalar left =
                row < 3 && innerRow < 3
              ? rotation(row, innerRow)
              : scalar(row == innerRow);

                for (direction innerColumn = 0; innerColumn < 4; ++innerColumn)
                {
                    const scalar right =
                        column < 3 && innerColumn < 3
                      ? rotation(column, innerColumn)
                      : scalar(column == innerColumn);

                    result(row, column) +=
                        left*value(innerRow, innerColumn)*right;
                }
            }
        }
    }

    return result;
}


//- Inverse-transform the spatial dimensions, leaving the fourth axis invariant
template<class Cmpt>
FOAM_DEVICE_I Tensor4<Cmpt> invTransform
(
    const tensor& rotation,
    const Tensor4<Cmpt>& value
)
{
    return transform(rotation.T(), value);
}

} // End namespace Foam

#include "Tensor4I.H"

#endif

// ************************************************************************* //