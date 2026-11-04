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

#include "IOstream.H"
#include "fe_coordinateSystem.H"
#include "fe_coordinateSystems.H"
#include "boundBox.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fe
{
    defineTypeNameAndDebug(coordinateSystem, 0);
    defineRunTimeSelectionTable(coordinateSystem, dictionary);
    defineRunTimeSelectionTable(coordinateSystem, origRotation);
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fe::coordinateSystem::coordinateSystem()
:
    name_(type()),
    note_(),
    origin_(point::zero),
    R_(),
    Rtr_(sphericalTensor::I)
{}


Foam::fe::coordinateSystem::coordinateSystem
(
    const word& name,
    const coordinateSystem& cs
)
:
    name_(name),
    note_(),
    origin_(cs.origin_),
    R_(cs.R_),
    Rtr_(cs.Rtr_)
{}


Foam::fe::coordinateSystem::coordinateSystem
(
    const word& name,
    const point& origin,
    const coordinateRotation& cr
)
:
    name_(name),
    note_(),
    origin_(origin),
    R_(cr),
    Rtr_(R_.T())
{}


Foam::fe::coordinateSystem::coordinateSystem
(
    const word& name,
    const point& origin,
    const vector& axis,
    const vector& dirn
)
:
    name_(name),
    note_(),
    origin_(origin),
    R_(axis, dirn),
    Rtr_(R_.T())
{}


Foam::fe::coordinateSystem::coordinateSystem
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    note_(),
    origin_(point::zero),
    R_(),
    Rtr_(sphericalTensor::I)
{
    operator=(dict);
}


Foam::fe::coordinateSystem::coordinateSystem
(
    const dictionary& dict
)
:
    name_(type()),
    note_(),
    origin_(point::zero),
    R_(),
    Rtr_(sphericalTensor::I)
{
    operator=(dict);
}


Foam::fe::coordinateSystem::coordinateSystem
(
    const dictionary& dict,
    const objectRegistry& obr
)
:
    name_(type()),
    note_(),
    origin_(point::zero),
    R_(),
    Rtr_(sphericalTensor::I)
{
    const entry* entryPtr = dict.findEntry(typeName_(), keyType::LITERAL);

    // a simple entry is a lookup into global coordinateSystems
    if (entryPtr && !entryPtr->isDict())
    {
        word csName;
        entryPtr->stream() >> csName;

        const coordinateSystems& csLst = coordinateSystems::New(obr);

        label csId = csLst.find(csName);
        if (debug)
        {
            Info<< "coordinateSystem::coordinateSystem"
                "(const dictionary&, const objectRegistry&):"
                << nl << "using global coordinate system: "
                << csName << "=" << csId << endl;
        }

        if (csId < 0)
        {
            FatalErrorIn
            (
                "coordinateSystem::coordinateSystem"
                "(const dictionary&, const objectRegistry&)"
            )   << "could not find coordinate system: " << csName << nl
                << "available coordinate systems: " << csLst.toc() << nl << nl
                << exit(FatalError);
        }

        // copy coordinateSystem, but assign the name as the typeName
        // to avoid strange things in writeDict()
        operator=(csLst[csId]);
        name_ = typeName_();
    }
    else
    {
        operator=(dict);
    }
}


Foam::fe::coordinateSystem::coordinateSystem(Istream& is)
:
    name_(is),
    note_(),
    origin_(point::zero),
    R_(),
    Rtr_(sphericalTensor::I)
{
    dictionary dict(is);
    operator=(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fe::coordinateSystem::~coordinateSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fe::coordinateSystem::spanInfo Foam::fe::coordinateSystem::spanLimited() const
{
    return spanInfo(Pair<bool>(false, false));
}


Foam::boundBox Foam::fe::coordinateSystem::spanBounds() const
{
    return boundBox::greatBox;
}


Foam::dictionary Foam::fe::coordinateSystem::dict(bool ignoreType) const
{
    dictionary dict;

    dict.add("name", name_);

    // only write type for derived types
    if (!ignoreType && type() != typeName_())
    {
        dict.add("type", type());
    }

    // The note entry is optional
    if (note_.size())
    {
        dict.add("note", note_);
    }

    dict.add("origin", origin_);
    dict.add("e1", e1());
    dict.add("e3", e3());

    return dict;
}


Foam::vector Foam::fe::coordinateSystem::localToGlobal
(
    const vector& local,
    bool translate
) const
{
    if (translate)
    {
        return (R_ & local) + origin_;
    }
    else
    {
        return (R_ & local);
    }
}


Foam::tmp<Foam::vectorField> Foam::fe::coordinateSystem::localToGlobal
(
    const vectorField& local,
    bool translate
) const
{
    if (translate)
    {
        return (R_ & local) + origin_;
    }
    else
    {
        return (R_ & local);
    }
}


Foam::vector Foam::fe::coordinateSystem::globalToLocal
(
    const vector& global,
    bool translate
) const
{
    if (translate)
    {
        return (Rtr_ & (global - origin_));
    }
    else
    {
        return (Rtr_ & global);
    }
}


Foam::tmp<Foam::vectorField> Foam::fe::coordinateSystem::globalToLocal
(
    const vectorField& global,
    bool translate
) const
{
    if (translate)
    {
        return (Rtr_ & (global - origin_));
    }
    else
    {
        return (Rtr_ & global);
    }
}


void Foam::fe::coordinateSystem::write(Ostream& os) const
{
    os  << type()
        << " origin: " << origin()
        << " e1: " << e1() << " e3: " << e3();
}


void Foam::fe::coordinateSystem::writeDict(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os  << indent << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    // only write type for derived types
    if (type() != typeName_())
    {
        os.writeKeyword("type") << type() << token::END_STATEMENT << nl;
    }

    // The note entry is optional
    if (note_.size())
    {
        os.writeKeyword("note") << note_ << token::END_STATEMENT << nl;
    }

    os.writeKeyword("name")   << name_    << token::END_STATEMENT << nl;
    os.writeKeyword("origin") << origin_  << token::END_STATEMENT << nl;
    os.writeKeyword("e1")     << e1()     << token::END_STATEMENT << nl;
    os.writeKeyword("e3")     << e3()     << token::END_STATEMENT << nl;

    if (subDict)
    {
        os << decrIndent << indent << token::END_BLOCK << endl;
    }
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::fe::coordinateSystem::operator=(const dictionary& rhs)
{
    if (debug)
    {
        Pout<< "coordinateSystem::operator=(const dictionary&) : "
            << "assign from " << rhs << endl;
    }

    // allow as embedded sub-dictionary "coordinateSystem"
    const dictionary& dict =
    (
        rhs.found(typeName_())
      ? rhs.subDict(typeName_())
      : rhs
    );

    // unspecified origin is (0 0 0)
    origin_ = point::zero;
    dict.readIfPresent("origin", origin_);

    // The note entry is optional
    note_.clear();
    rhs.readIfPresent("note", note_);

    // specify via coordinateRotation sub-dictionary
    if (dict.found("coordinateRotation"))
    {
        R_  = coordinateRotation::New(dict.subDict("coordinateRotation"))();
    }
    else
    {
        // let coordinateRotation constructor extract the axes specification
        R_ = coordinateRotation(dict);
    }

    Rtr_ = R_.T();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //
bool Foam::fe::operator!=(const Foam::fe::coordinateSystem& a, const Foam::fe::coordinateSystem& b)
{
    return (a.origin() != b.origin() || a.R() != b.R() || a.type() != b.type());
}

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::fe::operator<<(Ostream& os, const Foam::fe::coordinateSystem& cs)
{
    cs.write(os);
    os.check("Ostream& operator<<(Ostream&, const coordinateSystem&");
    return os;
}


// ************************************************************************* //
