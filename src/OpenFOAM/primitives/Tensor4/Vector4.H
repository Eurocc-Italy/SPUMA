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

\*---------------------------------------------------------------------------*/

#ifndef Foam_Vector4_H
#define Foam_Vector4_H

#include "contiguous.H"
#include "VectorSpace.H"

namespace Foam
{

template<class Cmpt>
class Vector4
:
    public VectorSpace<Vector4<Cmpt>, Cmpt, 4>
{
public:

    // Typedefs

        //- Equivalent type of labels used for valid component indexing
    typedef Vector4<label> labelType;


    // Member Constants

        //- Rank of Vector4 is 1
    static constexpr direction rank = 1;


    //- Component labeling enumeration
    enum components { X, Y, Z, K };


    // Generated Methods

        //- Default construct
    FOAM_DEVICE Vector4() = default;

        //- Copy construct
    FOAM_DEVICE Vector4(const Vector4&) = default;

        //- Copy assignment
    FOAM_DEVICE Vector4& operator=(const Vector4&) = default;


    // Constructors

        //- Construct initialized to zero
    FOAM_DEVICE inline Vector4(Foam::zero)
    :
        Vector4::vsType(Foam::zero{})
    {}

        //- Construct from a VectorSpace with four components
    template<class Cmpt2>
    FOAM_DEVICE inline explicit Vector4
    (
        const VectorSpace<Vector4<Cmpt2>, Cmpt2, 4>& vectorSpace
    )
    :
        Vector4::vsType(vectorSpace)
    {}

        //- Construct from four components
    FOAM_DEVICE inline Vector4
    (
        const Cmpt& vx,
        const Cmpt& vy,
        const Cmpt& vz,
        const Cmpt& vk
    )
    {
        x() = vx;
        y() = vy;
        z() = vz;
        k() = vk;
    }

        //- Construct from Istream
    inline explicit Vector4(Istream& is)
    :
        Vector4::vsType(is)
    {}


    // Member Functions

    // Component Access

        //- Access the vector x component
    FOAM_DEVICE const Cmpt& x() const noexcept { return this->v_[X]; }

        //- Access the vector y component
    FOAM_DEVICE const Cmpt& y() const noexcept { return this->v_[Y]; }

        //- Access the vector z component
    FOAM_DEVICE const Cmpt& z() const noexcept { return this->v_[Z]; }

        //- Access the vector k component
    FOAM_DEVICE const Cmpt& k() const noexcept { return this->v_[K]; }

        //- Access the vector x component
    FOAM_DEVICE Cmpt& x() noexcept { return this->v_[X]; }

        //- Access the vector y component
    FOAM_DEVICE Cmpt& y() noexcept { return this->v_[Y]; }

        //- Access the vector z component
    FOAM_DEVICE Cmpt& z() noexcept { return this->v_[Z]; }

        //- Access the vector k component
    FOAM_DEVICE Cmpt& k() noexcept { return this->v_[K]; }
};


// * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<Vector4<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous labels if component type is a contiguous label
template<class Cmpt>
struct is_contiguous_label<Vector4<Cmpt>> : is_contiguous_label<Cmpt> {};

//- Data are contiguous scalars if component type is a contiguous scalar
template<class Cmpt>
struct is_contiguous_scalar<Vector4<Cmpt>> : is_contiguous_scalar<Cmpt> {};

//- The transpose of a Vector4 is a Vector4
template<class Cmpt>
class typeOfTranspose<Cmpt, Vector4<Cmpt>>
{
public:
    typedef Vector4<Cmpt> type;
};

//- Use solveScalar components when solving for a Vector4
template<class Cmpt>
class typeOfSolve<Vector4<Cmpt>>
{
public:
    typedef Vector4<solveScalar> type;
};

} // End namespace Foam

#endif

// ************************************************************************* //