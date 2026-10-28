/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Application
    reductionSumTest

Group
    test/device

Description
    Test to check the reduction operation. 
    It is hard coded for an array size of 100000.

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "argList.H"
#include "Field.H"
#include "Random.H"
#include <tuple>
#include <chrono>
using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<typename T>
T testReductionSum(T* input, label size)
{
    T result = Zero;

    for (label i = 0; i < size; ++i)
    {
        result += input[i];
    }

    return result;
}

template <typename T>
struct type_tag
{
    using type = T;
};

template <typename F, typename... Ts>
void iterate_types(F &&func, std::tuple<Ts...> &&)
{
    (func(type_tag<Ts>{}), ...);
}

template<class T>
T randomValue(Random& rndGen)
{
    T value(Zero);
    for (direction cmpt = 0; cmpt < pTraits<T>::nComponents; ++cmpt)
    {
        setComponent(value, cmpt) =
            static_cast<typename pTraits<T>::cmptType>(rndGen.GaussNormal<scalar>());
    }
    return value;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Test for array reduction sum operation"
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "initDevice.H"
    #include "createMemoryPool.H"

    // Include edge cases for warp,block and device wide prefix sum
    List<label> sizes = {0,1,8,32,256,1024,8192,262144,1048576};

    using tuple_t =std::tuple<
        label,floatScalar,doubleScalar,
        Vector<label>,Vector<floatScalar>,Vector<doubleScalar>,
        Tensor<label>,Tensor<floatScalar>,Tensor<doubleScalar>
        >;
    tuple_t testType{};
    
    label Npassed = 0;
    const label numTypes = std::tuple_size<tuple_t>::value; // Number of types being tested
    const label Ntests = sizes.size()*numTypes;

    // Seed for rng
    auto now = std::chrono::high_resolution_clock::now();
    auto seed = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();

    auto testTypes = [&](auto typeTag)
    {
        using T = typename  decltype(typeTag)::type;
        Info << "Testing type: " << pTraits<T>::typeName << nl;

        for (label size : sizes)
        {
            Field<T> input(size);

            // Initialize input with random values
            Foam::Random rndGen(seed);
            for (size_t i = 0; i < size; ++i)
            {
                input[i] = randomValue<T>(rndGen);
            }

            // Compute sum
            const T computedSum = sum(input);

            // Compute expected prefix sum
            const T expectedSum = testReductionSum(input.begin(), size);

            // Compare output with expected
            bool passed = true;
            const scalar componentTolerance = 1e-4;
            const scalar tolerance =
                componentTolerance * sqrt(pTraits<T>::nComponents);

            if (mag(mag(expectedSum) - mag(computedSum))/(mag(expectedSum) + 1.e-10) > tolerance)
            {
                passed = false;
                Info << "\tMismatch in reduction sum result"
                     << ": expected " << expectedSum
                     << ", got " << computedSum << nl;
            }

            Info << "\treductionSum test, size: " << size << " on type: " <<
                pTraits<T>::typeName << " : " <<
                (passed ? "PASSED" : "FAILED") << nl;

            if (passed) ++Npassed;
        }
    };

    iterate_types(testTypes, std::move(testType));

    // Summary of all validation tests
    Info << "\n=== reductionSumTest execution summary: ===" << nl;
    Info << "\t testPassed: ["<< Npassed << " / " << Ntests << "]" << nl;

    return 0;
}

// ************************************************************************* //
