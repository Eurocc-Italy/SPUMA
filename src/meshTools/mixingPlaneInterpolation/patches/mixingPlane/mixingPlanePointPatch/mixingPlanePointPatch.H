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
    mixingPlanePointPatch

Description
    Dummy mixingPlane patch for post-processing.  No functionality built in!

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    mixingPlanePointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef mixingPlanePointPatch_H
#define mixingPlanePointPatch_H

#include "coupledFacePointPatch.H"
#include "mixingPlanePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class mixingPlanePointPatch Declaration
\*---------------------------------------------------------------------------*/

class mixingPlanePointPatch
:
    public coupledFacePointPatch
{
    // Private data

        //- Local reference cast into the cyclic patch
        const mixingPlanePolyPatch& mixingPlanePolyPatch_;


    // Private Member Functions

        //- Disallow default construct as copy
        mixingPlanePointPatch(const mixingPlanePointPatch&);

        //- Disallow default assignment
        void operator=(const mixingPlanePointPatch&);


        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&)
        {}

        //- Calculate mesh points
        virtual void calcGeometry(PstreamBuffers&)
        {}

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&)
        {}

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&)
        {}

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

    typedef pointBoundaryMesh BoundaryMesh;


    //- Runtime type information
    TypeName(mixingPlanePolyPatch::typeName_());


    // Constructors

        //- Construct from components
        mixingPlanePointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    // Destructor

        virtual ~mixingPlanePointPatch();


    // Member Functions

        // Access

            //- Return true because this patch is coupled
            virtual bool coupled() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
