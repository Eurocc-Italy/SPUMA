/*---------------------------------------------------------------------------*\\
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
    Foam::algebraicPairLduGAMGAgglomeration

Description
    Native serial algebraic-pair agglomeration for coupled LduMatrix GAMG.

SourceFiles
    algebraicPairLduGAMGAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef algebraicPairLduGAMGAgglomeration_H
#define algebraicPairLduGAMGAgglomeration_H

#include "LduGAMGAgglomeration.H"

namespace Foam
{

class algebraicPairLduGAMGAgglomeration
:
    public LduGAMGAgglomeration
{
    label mergeLevels_;
    bool forward_;

    tmp<labelField> pairAgglomerate
    (
        label& nCoarseCells,
        const lduAddressing& fineMatrixAddressing,
        const scalarField& faceWeights
    );

    void agglomerate(const scalarField& startFaceWeights);

public:

    static const word& typeName();

    algebraicPairLduGAMGAgglomeration
    (
        const lduMesh& mesh,
        const scalarField& faceWeights,
        const dictionary& controlDict
    );

    virtual const word& type() const;
};

} // End namespace Foam

#endif

// ************************************************************************* //