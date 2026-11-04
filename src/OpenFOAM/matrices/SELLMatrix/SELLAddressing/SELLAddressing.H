/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    Foam::SELLAddressing

Description
    To execute operations efficiently, SELL-C-sigma implements a structured 
    addressing layout:

    1) Local Padding: 
    Instead of padding every row to the global maximum row length 
    (like standard ELLPACK), SELL-C-sigma calculates the maximum row length 
    independently for each slice (of lenght C) and pads only up to that local maximum.

    2) Data Arrays:
    colsPtr_: stores the matching column indices mapped identically to the values array.   
    sliceBeginsPtr & slicesWidthPtr: auxiliary arrays that track the starting memory offset 
    and active width of each slice, allowing fast traversal during kernel execution.

SourceFiles
    SELLAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef SELLAddressing_H
#define SELLAddressing_H

#include "labelList.H"
#include "Tuple2.H"
#include "lduAddressing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SELLAddressing Declaration
\*---------------------------------------------------------------------------*/

class SELLAddressing
{
    // Private Data

        //- Number of equations/rows
        mutable label numRows_;

        //- Number of entries (after padding)
        mutable label numEntries_;

        //- Vector length of each slice
        const label C_;

        //- Number of slices
        mutable label numSlices_;


    //- Demand-driven data

        //- Map lower values into SELL values
        mutable std::unique_ptr<labelList> lowerToSELLPtr_;

        //- Map upper values into SELL values
        mutable std::unique_ptr<labelList> upperToSELLPtr_;

        //- Prefix-sum of elements per slice
        mutable std::unique_ptr<labelList> sliceBeginsPtr_;

        //- Padded width for each slice
        mutable std::unique_ptr<labelList> sliceWidthsPtr_;

        //- Packed SELL column indices, padded with INT32_MAX
        mutable std::unique_ptr<labelList> colsPtr_;


    // Private Member Functions

        //- Clear addressing
        void clearOut();


public:

    // Generated Methods

        //- No copy construct
        SELLAddressing(const SELLAddressing&) = delete;

        //- No copy assignment
        void operator=(const SELLAddressing&) = delete;


    // Constructors

        //- Default constructor
        SELLAddressing();

        //- Construct with number of equations/rows,
        // number of entries and lduAddressing
        SELLAddressing
        (
            const lduAddressing& addr
        );


    //- Destructor
    virtual ~SELLAddressing() 
    {
        clearOut();    
    };


    // Member Functions

        //- Calculate all SELL addressing arrays from an lduAddressing
        void calcAddr(const lduAddressing& addr) const;

        //- Return number of equations/rows
        label numRows() const noexcept
        {
            return numRows_;
        }

        //- Return number of entries
        label numEntries() const noexcept
        {
            return numEntries_;
        }

        //- Return the vector length
        label C() const noexcept
        {
            return C_;
        }

        //- Return the number of slices
        label numSlices() const noexcept
        {
            return numSlices_;
        }

        //- Return column indices
        const labelUList& cols() const
        {
            return *colsPtr_;
        }

        //- Return lowerToSELL
        const labelUList& lowerToSELL() const
        {
            return *lowerToSELLPtr_;
        }

        //- Return upperToSELL
        const labelUList& upperToSELL() const
        {
            return *upperToSELLPtr_;
        }

        //- Return slicePtr
        const labelUList& sliceBegins() const
        {
            return *sliceBeginsPtr_;
        }

        //- Return sliceWidthPtr
        const labelUList& sliceWidths() const
        {
            return *sliceWidthsPtr_;
        }

        //- Return info on the amount of padding
        void printInfo(const lduAddressing& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
