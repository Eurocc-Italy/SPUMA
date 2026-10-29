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

    Copyright (C) 2023 Martin Beaudoin

Class
    Foam::interfaceCheckFunctionObject

Description
    FunctionObject checking the fluxes across a pair of ggi, AMI or mixingPlane
    patches

Author
    Martin Beaudoin, 2023.

SourceFiles
    interfaceCheckFunctionObject.C

\*---------------------------------------------------------------------------*/


#ifndef interfaceCheckFunctionObject_H
#define interfaceCheckFunctionObject_H

#include "fvMeshFunctionObject.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class interfaceCheckFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class interfaceCheckFunctionObject
:
    public functionObject
{
    // Private data

        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

        //- Name of flux field
        word phiName_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        interfaceCheckFunctionObject(const interfaceCheckFunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const interfaceCheckFunctionObject&);


public:

    //- Runtime type information
    TypeName("interfaceCheck");


    // Constructors

        //- Construct from components
        interfaceCheckFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

    //- Destructor
    virtual ~interfaceCheckFunctionObject() = default;

    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute();

        //- Called when time was set at the end of the Time::operator++
        virtual bool timeSet();

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Write the registered objects
        virtual bool write()
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
