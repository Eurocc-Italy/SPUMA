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
    Foam::MemoryPoolManager

Description
    Manager class to handle arbitrary number of memory pools.

SourceFiles
    MemoryPoolManager.C

\*---------------------------------------------------------------------------*/

#ifndef Memory_Pool_Manager_H
#define Memory_Pool_Manager_H

#include <vector>
#include "MemoryPool.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

/*---------------------------------------------------------------------------*\
                           Class MemoryPoolManager Declaration
\*---------------------------------------------------------------------------*/

class MemoryPoolManager
{
private:

    // Private Member Data

        //- List of memory pool objects.
	static std::vector<MemoryPool*> pools_;


public:


    // Static Data
        
        TypeName("MemoryPoolManager");


private:

    // Constructors

        //- Constructor
        MemoryPoolManager
        (
            const word arch,
            const word defaultPoolType,
            const scalar defaultPoolSize
        );

        //- No copy constructor
        MemoryPoolManager(MemoryPoolManager &other) = delete;


public:

    //- Destructor
    ~MemoryPoolManager() = default;


    // Static Member Functions

        //- Static access method to build the memory pool manager
        static MemoryPoolManager& New
        (
            const word arch,
            const word defaultPoolType,
            const scalar defaultPoolSize
        );

        //- Static access method to get the default MemoryPool instance
        static MemoryPool* getMemoryPoolInstance();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

