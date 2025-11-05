/*---------------------------------------------------------------------------*\
  *      .  *_______ * ______ .  __ *  __ * ___ .___    .  ___ .   *  .     *
    *  .    /       | |   _  \  |  |  |  | |   \/   | *   /   \ *   .    *   .
 *    .  * .\   (---*.|  |_)  |.|  |  |  |*|  \  /  |. * /  *  \  .  *     *
 =^^=^^==^^^=\   \^=^=|   ___/=^|  |^=|  |=|  |\/|  |^^=/  /=\  \^=^=^^===^^^=
 0  o  O  o---)   \ 0 |  |   0  |  o--o  |o|  |  |  | o/  _____  \ 0   o  O
     0    |_______/   |__| o   o \______/  |__| 0|__| /__/  o  \__\   o
  O   o  o        0  o      0   O        o    o       O  o     0   o    0  o
-------------------------------------------------------------------------------
    Copyright (C) 2025 Cineca
-------------------------------------------------------------------------------
License
    This file is part of SPUMA.

    SPUMA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SPUMA is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with SPUMA.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::MemoryPool

Description
    Frontend interface to memory pool implementations.
    Memory pools are algorithms that rely on pools to
    manage dynamic memory allocation in a efficient way.
    Memory pools remove the memory fragmentation issues and allow memory
    allocation (and deallocation) with constant execution time

SourceFiles
    MemoryPool.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Memory_Pool_H
#define Foam_Memory_Pool_H

#include "label.H"
#include "word.H"
#include "scalar.H"
#include "typeInfo.H"
#include <map>

//- macro to raise poolptr errors
#define raisePoolValidError(poolPtr)\
    FatalErrorInFunction\
        << "MEMPOOL: pointer " << std::to_string(reinterpret_cast<uint64_t>(poolPtr))\
        << " is not valid" << abort(FatalError);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class MemoryPool;

/*---------------------------------------------------------------------------*\
                           Class MemoryPool Declaration
\*---------------------------------------------------------------------------*/

class MemoryPool
{
protected:

    // Private data

        //- Size of the pool
        uint64_t size_;
    
        //- The size of the allocated memory
        uint64_t allocatedSize_;
    
        //- The size of the unallocated memory
        uint64_t unallocatedSize_;
    
        //- Max occupancy reached during the simulation
        uint64_t maxOccupancy_;

        //- Useful typedefs
        typedef std::map<char*, uint64_t> blockList;
        typedef std::pair<char*, uint64_t> blockPair;

        //- List of unused memory blocks 
        // (represented in memory as a map)
        blockList unusedBlockList_;
   
        //- List of used memory blocks
        // (represented in memory as a map)
        blockList usedBlockList_;

public:
    // Static Data
    
        // Typename
        TypeName("MemoryPool");


    // Constructors
    
        //- Default constructor
        MemoryPool();

        //- Memory pool constructor from dictionary
        MemoryPool(const dictionary& dict);

        //- No copy constructor
        MemoryPool(MemoryPool &other) = delete;


    //- Destructor
    virtual ~MemoryPool();


    // Static Member Functions

        //- Create a new memory pool of type @type
        // and size @size
        static MemoryPool* New
        (
            const word arch,
            const word type,
            const uint64_t size
        );


    // Member Functions

        //- Allocate a new block of memory of size @size 
        virtual void* allocate(uint64_t size) = 0;

        //- Free the block of memory associated to @ptr
        virtual void free(void* ptr) = 0;

        //- Check if pool pointer @poolPtr is valid
        virtual inline bool isValid(void* poolPtr) = 0;

        //- Check if @ptr is inside a block
        virtual inline bool isInBlockRange(void* ptr) = 0;

        //- Return the size of the array associated to @poolPtr 
        // (in bytes if type is not specified)
        virtual uint64_t arraySizeInBytes(void* poolPtr) = 0;

        //- Copy data in the memory pool
        virtual void copyIn
        (
            void* poolPtr,
            void* ptr,
            uint64_t nElementsInBytes
        ) = 0;

        //- Copy data from the memory pool
        virtual void copyOut
        (
            void* poolPtr,
            void* ptr,
            uint64_t nElementsInBytes
        ) = 0;

        //- Set the memory to the scalar value @value 
        // for each byte of specified memory
        virtual void memSet
        (
            void* poolPtr,
            const void* value,
            size_t sizeOfValue,
            uint64_t nElementsInBytes
        ) = 0;

        //- Set the memory to a value of one 
        // for each byte of specified memory
        virtual void memSetScalarOne
        (
            void* poolPtr,
            uint64_t nElementsInBytes
        ) = 0;

        //- Set the memory to the integer value @value 
        // using built-in memset function
        virtual void memSet
        (
            void* poolPtr,
            const int value,
            uint64_t nElementsInBytes
        ) = 0;

        //- Copy from one block of the memory pool to another
        virtual void memCopy
        (
            void* tgtPtr,
            void* srcPtr,
            uint64_t nElementsInBytes
        ) = 0;

    // Pool status

        //- Return the pool size
        virtual uint64_t size() const = 0;

        //- Return the allocated size
        virtual uint64_t allocatedSize() const = 0;

        //- Return the unallocated size
        virtual uint64_t unallocatedSize() const = 0;

        //- Return the occupancy rate (in percent)
        virtual scalar occupancy() const
        {
            return 100 * scalar(this->allocatedSize()) / scalar(this->size());
        }

        //- Return the max occupancy reached
        // during the simulation
        virtual uint64_t maxOccupancy() const
        {
            return maxOccupancy_;
        }

        //- Print the list of allocated blocks
        virtual void showAllocated(bool relative = false) = 0;

        //- Print the list of unallocated blocks
        virtual void showUnallocated(bool relative = false) = 0;


    // Member Operators

        //- Remove the copy assignement operator
        virtual void operator=(const MemoryPool &) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
