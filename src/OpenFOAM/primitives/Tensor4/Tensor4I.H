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

template<class Cmpt>
FOAM_DEVICE inline Foam::Tensor4<Cmpt>::Tensor4(Foam::zero)
:
    Tensor4::msType(Foam::zero{})
{}

template<class Cmpt>
template<class Cmpt2>
FOAM_DEVICE inline Foam::Tensor4<Cmpt>::Tensor4
(
    const MatrixSpace<Tensor4<Cmpt2>, Cmpt2, 4, 4>& matrix
)
:
    Tensor4::msType(matrix)
{}

template<class Cmpt>
template<class Cmpt2>
FOAM_DEVICE inline Foam::Tensor4<Cmpt>::Tensor4
(
    const VectorSpace<Tensor4<Cmpt2>, Cmpt2, 16>& vectorSpace
)
:
    Tensor4::msType(vectorSpace)
{}

template<class Cmpt>
FOAM_DEVICE inline Foam::Tensor4<Cmpt>::Tensor4
(
    const Cmpt& txx, const Cmpt& txy, const Cmpt& txz, const Cmpt& txk,
    const Cmpt& tyx, const Cmpt& tyy, const Cmpt& tyz, const Cmpt& tyk,
    const Cmpt& tzx, const Cmpt& tzy, const Cmpt& tzz, const Cmpt& tzk,
    const Cmpt& tkx, const Cmpt& tky, const Cmpt& tkz, const Cmpt& tkk
)
{
    this->v_[XX] = txx; this->v_[XY] = txy;
    this->v_[XZ] = txz; this->v_[XK] = txk;
    this->v_[YX] = tyx; this->v_[YY] = tyy;
    this->v_[YZ] = tyz; this->v_[YK] = tyk;
    this->v_[ZX] = tzx; this->v_[ZY] = tzy;
    this->v_[ZZ] = tzz; this->v_[ZK] = tzk;
    this->v_[KX] = tkx; this->v_[KY] = tky;
    this->v_[KZ] = tkz; this->v_[KK] = tkk;
}

template<class Cmpt>
inline Foam::Tensor4<Cmpt>::Tensor4(Istream& is)
:
    Tensor4::msType(is)
{}

template<class Cmpt>
template<class Cmpt2>
FOAM_DEVICE inline void Foam::Tensor4<Cmpt>::operator=
(
    const VectorSpace<Tensor4<Cmpt2>, Cmpt2, 16>& vectorSpace
)
{
    VectorSpace<Tensor4<Cmpt>, Cmpt, 16>::operator=(vectorSpace);
}

template<class Cmpt>
FOAM_DEVICE inline Foam::Vector4<Cmpt> Foam::solve
(
    const Foam::Tensor4<Cmpt>& matrix,
    const Foam::Vector4<Cmpt>& rhs
)
{
    Cmpt augmented[4][5];
    scalar scale = 0;

    for (label row = 0; row < 4; ++row)
    {
        for (label column = 0; column < 4; ++column)
        {
            augmented[row][column] = matrix(row, column);
            scale = max(scale, mag(augmented[row][column]));
        }
    }

    for (label row = 0; row < 4; ++row)
    {
        augmented[row][4] = rhs[row];
    }

    const scalar threshold = ROOTVSMALL*scale;

    for (label column = 0; column < 4; ++column)
    {
        label pivotRow = column;
        scalar pivotMagnitude = mag(augmented[pivotRow][column]);

        for (label row = column + 1; row < 4; ++row)
        {
            const scalar candidateMagnitude = mag(augmented[row][column]);

            if (candidateMagnitude > pivotMagnitude)
            {
                pivotRow = row;
                pivotMagnitude = candidateMagnitude;
            }
        }

        if (pivotMagnitude <= threshold)
        {
            return Vector4<Cmpt>(Zero);
        }

        if (pivotRow != column)
        {
            for (label entry = column; entry < 5; ++entry)
            {
                const Cmpt value = augmented[column][entry];
                augmented[column][entry] = augmented[pivotRow][entry];
                augmented[pivotRow][entry] = value;
            }
        }

        const Cmpt pivot = augmented[column][column];

        for (label entry = column; entry < 5; ++entry)
        {
            augmented[column][entry] /= pivot;
        }

        for (label row = 0; row < 4; ++row)
        {
            if (row != column)
            {
                const Cmpt factor = augmented[row][column];

                for (label entry = column; entry < 5; ++entry)
                {
                    augmented[row][entry] -= factor*augmented[column][entry];
                }
            }
        }
    }

    return Vector4<Cmpt>
    (
        augmented[0][4],
        augmented[1][4],
        augmented[2][4],
        augmented[3][4]
    );
}


template<class Cmpt>
FOAM_DEVICE inline Foam::Tensor4<Cmpt> Foam::inv
(
    const Foam::Tensor4<Cmpt>& matrix
)
{
    Tensor4<Cmpt> result(Zero);

    for (direction column = 0; column < 4; ++column)
    {
        Vector4<Cmpt> unit(Zero);
        unit[column] = Cmpt(1);

        const Vector4<Cmpt> inverseColumn = solve(matrix, unit);

        for (direction row = 0; row < 4; ++row)
        {
            result(row, column) = inverseColumn[row];
        }
    }

    return result;
}

// ************************************************************************* //