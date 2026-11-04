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

\*---------------------------------------------------------------------------*/

#include "fe_coordinateSystems.H"
#include "IOPtrList.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fe
{
    defineTypeNameAndDebug(coordinateSystems, 0);
}
defineTemplateTypeNameAndDebug(IOPtrList<fe::coordinateSystem>, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fe::coordinateSystems::coordinateSystems(const IOobject& io)
:
    IOPtrList<coordinateSystem>(io)
{}


Foam::fe::coordinateSystems::coordinateSystems
(
    const IOobject& io,
    const PtrList<coordinateSystem>& lst
)
:
    IOPtrList<coordinateSystem>(io, lst)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// Read construct from registry, or return previously registered
const Foam::fe::coordinateSystems& Foam::fe::coordinateSystems::New
(
    const objectRegistry& obr
)
{
    if (obr.foundObject<coordinateSystems>(typeName))
    {
        return obr.lookupObject<coordinateSystems>(typeName);
    }
    else
    {
        return obr.store
        (
            new coordinateSystems
            (
                IOobject
                (
                    typeName,
                    "constant",
                    obr,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::fe::coordinateSystems::find(const word& keyword) const
{
    forAll(*this, i)
    {
        if (keyword == operator[](i).name())
        {
            return i;
        }
    }

    return -1;
}


bool Foam::fe::coordinateSystems::found(const word& keyword) const
{
    return find(keyword) >= 0;
}


Foam::wordList Foam::fe::coordinateSystems::toc() const
{
    wordList keywords(size());

    forAll(*this, i)
    {
        keywords[i] = operator[](i).name();
    }

    return keywords;
}


bool Foam::fe::coordinateSystems::writeData(Ostream& os, bool subDict) const
{
    os << nl << size() << nl << token::BEGIN_LIST;

    forAll(*this, i)
    {
        os << nl;
        operator[](i).writeDict(os, subDict);
    }

    os << token::END_LIST << nl;

    return os.good();
}


// ************************************************************************* //
