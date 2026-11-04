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

#include "fe_axisCoordinateRotation.H"

#include "Switch.H"
#include "mathematicalConstants.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace fe
{
    defineTypeNameAndDebug(axisCoordinateRotation, 0);
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        axisCoordinateRotation,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fe::axisCoordinateRotation::calcTransform
(
    const scalar phiAngle,
    const scalar thetaAngle,
    const scalar psiAngle,
    const bool inDegrees
)
{
    scalar phi   = phiAngle;
    scalar theta = thetaAngle;
    scalar psi   = psiAngle;

    if (inDegrees)
    {
        phi   *= constant::mathematical::pi/180.0;
        theta *= constant::mathematical::pi/180.0;
        psi   *= constant::mathematical::pi/180.0;
    }

    tensor::operator=
    (
        tensor
        (
            cos(theta)*cos(psi),
            sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi),
            cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi),

            cos(theta)*sin(psi),
            sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi),
            cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi),

            -sin(theta),
            sin(phi)*cos(theta),
            cos(phi)*cos(theta)
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fe::axisCoordinateRotation::axisCoordinateRotation
(
    const scalar phiAngle,
    const scalar thetaAngle,
    const scalar psiAngle,
    const bool inDegrees
)
:
    coordinateRotation()
{
    calcTransform(phiAngle, thetaAngle, psiAngle, inDegrees);
}


Foam::fe::axisCoordinateRotation::axisCoordinateRotation
(
    const dictionary& dict
)
:
    coordinateRotation()
{
    scalar phi = readScalar(dict.lookup("phi"));
    scalar theta = readScalar(dict.lookup("theta"));
    scalar psi = readScalar(dict.lookup("psi"));

    calcTransform
    (
        phi,
        theta,
        psi,
        dict.lookupOrDefault<Switch>("degrees", true)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
