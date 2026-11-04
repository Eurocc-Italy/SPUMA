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

Description
    Shared template name for MixingPlane interpolation

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "mixingPlaneInterpolationTemplate.H"
#include "mixingPlaneInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MixingPlaneInterpolationName, 0);

    defineTemplateTypeNameAndDebug(mixingPlaneZoneInterpolation, false);
    defineTemplateTypeNameAndDebug(ggiZoneInterpolation, false);
}

const Foam::Enum<Foam::MixingPlaneInterpolationName::discretizationMethod>
Foam::MixingPlaneInterpolationName::discretizationMethodNames_
({
    { discretizationMethod::MASTER_PATCH, "masterPatch" }, // Deprecated syntax
    { discretizationMethod::SLAVE_PATCH,  "slavePatch" },  // Deprecated syntax
    { discretizationMethod::BOTH_PATCHES, "bothPatches" }, // Deprecated syntax

    { discretizationMethod::FOLLOW_PATCH_DISCRETIZATION, "followPatchDiscretization" },
    { discretizationMethod::LINEAR_ALONG_STACK_AXIS,     "linearAlongStackAxis" },
    { discretizationMethod::POLYNOMIAL,                  "polynomial" },
    { discretizationMethod::USER_DEFINED,                "userDefinedProfile" },
    { discretizationMethod::DISCRETIZATION_UNKNOWN,      "unknown" }
});

const Foam::Enum<Foam::MixingPlaneInterpolationName::followPatchDiscretizationOption>
Foam::MixingPlaneInterpolationName::followPatchDiscretizationOptionNames_
({
    { followPatchDiscretizationOption::FOLLOW_MASTER_PATCH,    "masterPatch" },
    { followPatchDiscretizationOption::FOLLOW_NEIGHBOUR_PATCH, "neighbourPatch" },
    { followPatchDiscretizationOption::FOLLOW_BOTH_PATCHES,    "bothPatches" },
    { followPatchDiscretizationOption::FOLLOW_UNKNOWN,         "unknown" }
});

const Foam::Enum<Foam::MixingPlaneInterpolationName::userDefinedProfilecoordinateSystemType>
Foam::MixingPlaneInterpolationName::userDefinedProfilecoordinateSystemTypeNames_
({
    { userDefinedProfilecoordinateSystemType::USER_DEFINED_PROFILE_COORDS_LOCAL,  "local" },
    { userDefinedProfilecoordinateSystemType::USER_DEFINED_PROFILE_COORDS_GLOBAL, "global" },
    { userDefinedProfilecoordinateSystemType::USERDEFINEDPROFILECOORDINATESYSTEMTYPE_UNKNOWN, "unknown" }
});

const Foam::Enum<Foam::MixingPlaneInterpolationName::sweepAxis>
Foam::MixingPlaneInterpolationName::sweepAxisNames_
({
    { sweepAxis::SWEEP_X,       "X" },
    { sweepAxis::SWEEP_Y,       "Y" },
    { sweepAxis::SWEEP_Z,       "Z" },
    { sweepAxis::SWEEP_R,       "R" },
    { sweepAxis::SWEEP_THETA,   "Theta" },
    { sweepAxis::SWEEP_UNKNOWN, "Unknown" }
});

const Foam::Enum<Foam::MixingPlaneInterpolationName::stackAxis>
Foam::MixingPlaneInterpolationName::stackAxisNames_
({
    { stackAxis::STACK_X,       "X" },
    { stackAxis::STACK_Y,       "Y" },
    { stackAxis::STACK_Z,       "Z" },
    { stackAxis::STACK_R,       "R" },
    { stackAxis::STACK_THETA,   "Theta" },
    { stackAxis::STACK_UNKNOWN, "Unknown" }
});

const Foam::Enum<Foam::MixingPlaneInterpolationName::mixingType>
Foam::MixingPlaneInterpolationName::mixingTypeNames_
({
    { mixingType::AREA_AVERAGING,   "areaAveraging" },
    { mixingType::FLUX_AVERAGING,   "fluxAveraging" },
    { mixingType::UNIFORM_VALUE,    "uniformValue" },
    { mixingType::UNIFORM_GRADIENT, "uniformGradient" },
    { mixingType::ZERO_GRADIENT,    "zeroGradient" },
    { mixingType::MIXING_UNKNOWN,   "Unknown" }
});

const Foam::Enum<Foam::MixingPlaneInterpolationName::wfEvalMethod>
Foam::MixingPlaneInterpolationName::wfEvalMethodNames_
({
    { wfEvalMethod::GGI,                             "ggi" },
    { wfEvalMethod::AMI,                             "ami" },
    { wfEvalMethod::WFEVALMETHOD_UNKNOWN, "Unknown" }
});


// ************************************************************************* //
