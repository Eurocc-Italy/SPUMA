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

Class
    Foam::coordinateSystems

Description
    Provides a centralized coordinateSystem collection.

Note
    Mixing normal constructors and the coordinateSystems::New constructor
    may yield unexpected results.

SourceFiles
    coordinateSystems.C

\*---------------------------------------------------------------------------*/
#ifndef coordinateSystems_H
#define coordinateSystems_H

#include "fe_coordinateSystem.H"
#include "IOPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fe
{
/*---------------------------------------------------------------------------*\
                     Class coordinateSystems Declaration
\*---------------------------------------------------------------------------*/

class coordinateSystems
:
    public IOPtrList<fe::coordinateSystem>
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        coordinateSystems(const coordinateSystems&);

        //- Disallow default bitwise assignment
        void operator=(const coordinateSystems&);

public:

    //- Runtime type information
    TypeName("coordinateSystems");

    // Constructors

        //- Read construct from IOobject
        explicit coordinateSystems(const IOobject&);

        //- Construct from IOobject and a PtrList
        coordinateSystems
        (
            const IOobject&,
            const PtrList<coordinateSystem>&
        );


    // Selectors

        //- Return previously registered or read construct from "constant"
        static const coordinateSystems& New(const objectRegistry&);

    // Member Functions

        //- Find and return index for a given keyword, returns -1 if not found
        label find(const word& key) const;

        //- Search for given keyword
        bool found(const word& keyword) const;

        //- Return the table of contents (list of all keywords)
        wordList toc() const;

        using IOPtrList<fe::coordinateSystem>::writeData;

        //- write data
        bool writeData(Ostream&, bool subDict=true) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fe
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
