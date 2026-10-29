/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
    Copyright (C) 2025 CINECA (Ettore Fadiga)
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
    Foam::functionObjects::writeTurboFields

Group

Description
    It has similar functionality as the main time database through the
    \c writeControl setting:
      - timeStep
      - writeTime
      - adjustableRunTime
      - runTime
      - clockTime
      - cpuTime

Usage
    Example of function object specification:
    \verbatim
    writeTurboFields1
    {
        type        writeTurboFields;
        libs        ("libTurboPerformance.so");
        axialDir    (x y z);
        omega       1;  // [rad/s]
        origin      (x y z);
        ...
        fields      (Ucyl URelcyl UMer URelMer Utr h0 p0 Rothalpy);
    }
    \endverbatim

See also
    Foam::functionObject
    Foam::functionObjects::timeControl

SourceFiles
    writeTurboFields.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_writeTurboFields_H
#define Foam_functionObjects_writeTurboFields_H

#include "functionObject.H"
#include "wordRes.H"
#include "Enum.H"
#include "IOMRFZoneList.H"
#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                Class functionObjects::writeTurboFields Declaration
\*---------------------------------------------------------------------------*/

class writeTurboFields
:
    public fvMeshFunctionObject
{

private:

    // Private Data

        //- Reference to registry
        const objectRegistry& obr_;

        //- Names of objects to control
        wordRes objectNames_;

        //- Name of pressure field
        word pName_;

        //- Name of velocity field
        word UName_;

        //- Name of temperature field
        word TName_;

        //- Axial direction
        vector axialDirection_;

        //- Origin
        vector origin_;

        //- Rot speed [rad/s]
        scalar omega_;

        // MRF to calculate relative velocity
        IOMRFZoneList MRF_;

    // Private Member Functions

        //- No copy construct
        writeTurboFields(const writeTurboFields&) = delete;

        //- No copy assignment
        void operator=(const writeTurboFields&) = delete;


public:

    //- Runtime type information
    TypeName("writeTurboFields");


    // Constructors

        //- Construct from Time and dictionary
        writeTurboFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~writeTurboFields() = default;


    // Member Functions

        //- Read the writeTurboFields data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Write the registered objects
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
