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
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::Function1Types::twg_CSV

Description
    Templated CSV function.
    Modified from the standard FOAM::CSV class.

    The reference column is always a scalar, e.g. time.

    Usage:
    \verbatim
        <entryName> csvFile;
        <entryName>Coeffs
        {
            nHeaderLine         4;          // number of header lines (optional)
            refColumn           0;          // reference column index (optional)
            componentColumns    (1 2 3);    // component column indices
            separator           ",";        // optional (defaults to ",")
            commentChar         "#";        // optional (defaults to "#")
            mergeSeparators     false;      // merge multiple separators (defaults to false)
            file                "fileXYZ";  // name of csv data file
            outOfBounds         clamp;      // optional out-of-bounds handling
            interpolationScheme linear;     // optional interpolation scheme
        }
    \endverbatim

SourceFiles
    twg_CSV.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Function1Types_twg_CSV_H
#define Foam_Function1Types_twg_CSV_H

#include "Function1.H"
#include "TableBase.H"
#include "Tuple2.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class twg_CSV Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class twg_CSV
:
    public TableBase<Type>
{
    // Private Data

        //- Number header lines
        const label nHeaderLine_;

        //- Column of the time (optional)
        const label refColumn_;

        //- Labels of the components
        const labelList componentColumns_;

        //- Separator character
        const char separator_;

        //- Comment character
        mutable char comment_;

        //- Merge separators flag, e.g. ',,,' becomes ','
        const bool mergeSeparators_;

        //- File name for csv table
        fileName fName_;


    // Private Member Functions

        //- Get component columns entry
        static labelList getComponentColumns
        (
            const word& name,
            const dictionary& dict
        );

        //- Read csv data table
        void read();

        //- Read component values from the split string
        Type readValue(const UList<string>& strings) const;

public:

    // Generated Methods

        //- No copy assignment
        void operator=(const twg_CSV<Type>&) = delete;

    //- Declare type-name, virtual type (with debug switch)
    ClassName("twg_csvFile");

    // Constructors

        //- Construct from entry name, dictionary and optional registry
        twg_CSV
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr,
            const fileName& fName = fileName::null
        );

    //- Copy construct
        explicit twg_CSV(const twg_CSV<Type>& csv);

    //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const
        {
            return tmp<Function1<Type>>(new twg_CSV<Type>(*this));
        }

    //- Destructor
    virtual ~twg_CSV() = default;


    // Member Functions

        //- Return const access to the file name
        virtual const fileName& fName() const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;

        //- Write coefficient entries in dictionary format
        virtual void writeEntries(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "twg_CSV.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
