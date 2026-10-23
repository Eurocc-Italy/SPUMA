/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::fvMeshSELLAddressing

Description
    SELL-C-sigma representation of the mesh addressing.

SourceFiles
    fvMeshSELLAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvMeshSELLAddressing_H
#define Foam_fvMeshSELLAddressing_H

#include "SELLAddressing.H"
#include "labelList.H"
#include "UList.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                        Class fvMeshSELLddressing Declaration
\*---------------------------------------------------------------------------*/

class fvMeshSELLAddressing
: 
    public SELLAddressing
{
    // Private Member Functions

        //- No copy construct
        fvMeshSELLAddressing(const fvMeshSELLAddressing&) = delete;

        //- No copy assignment
        void operator=(const fvMeshSELLAddressing&) = delete;

public:

    // Constructors

        //- Construct from mesh
        fvMeshSELLAddressing(const fvMesh& mesh);


    //- Destructor
        ~fvMeshSELLAddressing() = default;
 

    // Member Functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
