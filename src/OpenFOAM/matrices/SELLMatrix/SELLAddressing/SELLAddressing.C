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
    Copyright (C) 2026 Gabriele Ceccolini - Unibo
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

\*---------------------------------------------------------------------------*/

#include "SELLAddressing.H"
#include "scalarField.H"
#include "deviceInits.H"
#include "registerSwitch.H"
#include "debug.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#if defined(have_cuda) || defined(have_hip)
// Not modifiable by the user
Foam::label Foam::SELLAddressing::C(Foam::foamDeviceInit::getNumberOfThreadsPerBlock());
#else
// User-defined value
Foam::label Foam::SELLAddressing::C
(
    Foam::debug::optimisationSwitch("sellC", 1)
);
registerOptSwitch
(
    "sellC",
    Foam::label,
    Foam::SELLAddressing::C
);
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SELLAddressing::SELLAddressing()
:
    numRows_(0),
    numEntries_(0),
    numSlices_(0)
{}


Foam::SELLAddressing::SELLAddressing
(
    const lduAddressing& addr
)
:
    numRows_(0),
    numEntries_(0),
    numSlices_(0)
{
    calcAddr(addr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SELLAddressing::calcAddr(const lduAddressing& addr) const
{
    if (numRows_ > 0)
        return;

    const label sigma = 1;  // no sorting
    const label nCells = addr.size();
    const label nInternalFaces = addr.upperAddr().size();
    numSlices_ = (nCells + C - 1) / C;
    numRows_ = nCells;

    const labelUList& lowerAddr = addr.lowerAddr();
    const labelUList& upperAddr = addr.upperAddr();

    // Each face contributes one off-diagonal in each of its two rows,
    // no diagonal contribution
    labelList nnz(nCells, 0);
    
    forAll(lowerAddr, facei)
    {
        ++nnz[lowerAddr[facei]];
        ++nnz[upperAddr[facei]];
    }

    // Allocate sliceWidthsPtr and fill in
    sliceWidthsPtr_ = std::make_unique<labelList>
    (
        numSlices_,
        Foam::zero{},
        poolSwitch(1)
    );
    auto& sliceWidths = *sliceWidthsPtr_;

    for (label s = 0; s < numSlices_; ++s) 
    {
        const label rowStart = s * C;
        const label rowEnd = (rowStart + C > nCells) ? nCells : rowStart + C;
        label maxWidth = 0;
        
        for (label rowi = rowStart; rowi < rowEnd; ++rowi) 
        {
            const label nnzi = nnz[rowi];
            if (nnzi > maxWidth) 
            {
                maxWidth = nnz[rowi];
            }
        }
        sliceWidths[s] = maxWidth;
    }

    // Allocate sliceBegins and fill in
    sliceBeginsPtr_ = std::make_unique<labelList>
    (
        numSlices_ + 1,
        Foam::zero{},
        poolSwitch(1)
    );
    auto& sliceBegins = *sliceBeginsPtr_;

    sliceBegins[0] = 0;
    for (label s = 0; s < numSlices_; ++s)
    {
        sliceBegins[s + 1] = sliceBegins[s] + sliceWidths[s] * C;
    }
    numEntries_ = sliceBegins[numSlices_];

    // Allocate lowetToSELLPtr
    lowerToSELLPtr_ = std::make_unique<labelList>
    (
        nInternalFaces,
        Foam::zero{},
        poolSwitch(1)
    );
    auto& lowerToSELL = *lowerToSELLPtr_;

    // Allocate upperToSELLPtr
    upperToSELLPtr_ = std::make_unique<labelList>
    (
        nInternalFaces,
        Foam::zero{},
        poolSwitch(1)
    );
    auto& upperToSELL = *upperToSELLPtr_;

    // Allocate colsPtr and initialization
    colsPtr_ = std::make_unique<labelList>
    (
        numEntries_,
        Foam::zero{},
        poolSwitch(1)
    );
    auto& cols = *colsPtr_;

    // Initialize with padding value
    for (label i = 0; i < numEntries_; ++i) 
    {
        cols[i] = INT32_MAX;
    }

    labelList currentSlot(nCells, 0);

    // Fill in cols and conversion maps
    for (label facei=0; facei < nInternalFaces; ++facei) 
    {
        label i = lowerAddr[facei];
        label j = upperAddr[facei];

        // Entry (i, j)
        {
            const label slicei = i / C;                                 // Slice index
            const label sliceRowi = i % C;                              // Local row index in slice
            const label sliceBegin = sliceBegins[slicei];
            const label sloti = currentSlot[i]++;
            const label slotRowi = sliceBegin + sloti * C + sliceRowi;  // Column-major layout within slice

            cols[slotRowi] = j;
            upperToSELL[facei] = slotRowi; 
        }

        // Symmetric entry (j, i)
        {
            const label slicei = j / C;                                // Slice index
            const label sliceRowi = j % C;                             // Local row index in slice
            const label sliceBegin = sliceBegins[slicei];
            const label sloti = currentSlot[j]++;
            const label slotRowi = sliceBegin + sloti * C + sliceRowi; // Column-major layout within slice

            cols[slotRowi] = i;
            lowerToSELL[facei] = slotRowi;
        }
    }

    printInfo(addr);
}


void Foam::SELLAddressing::clearOut()
{
    numRows_ = 0;
    numEntries_ = 0;
    numSlices_ = 0;

    lowerToSELLPtr_.reset(nullptr);
    upperToSELLPtr_.reset(nullptr);
    sliceBeginsPtr_.reset(nullptr);
    sliceWidthsPtr_.reset(nullptr);
    colsPtr_.reset(nullptr);
}


void Foam::SELLAddressing::printInfo(const lduAddressing& addr) const
{
    const labelUList& lowerAddr = addr.lowerAddr();
    const label nCells = numRows();
    auto& sliceWidths = *sliceWidthsPtr_;

    label origNnz = 2 * lowerAddr.size();
    label newNnz = 0;

    for (label s = 0; s < numSlices_; ++s)
    {
        const label rowStart = s * C;
        const label rowEnd = (rowStart + C > nCells) ? nCells : rowStart + C;
        newNnz += sliceWidths[s] * (rowEnd - rowStart);
    }

    Info << "SELLAddressing info:" << nl
         << "    original nnz: " << origNnz << nl
         << "    current nnz (after padding): " << newNnz << nl
         << "    paddingFactor: " << 100.* scalar(newNnz - origNnz)/scalar(origNnz) << nl
         << "    chunk occupancy (beta): " << scalar(origNnz) / scalar(newNnz) << nl;
}


// ************************************************************************* //
