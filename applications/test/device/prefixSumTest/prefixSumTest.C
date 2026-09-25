/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    prefixSumTest

Group
    test/device

Description
    Test to check the prefix sum operation.

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "argList.H"
#include "Field.H"
#include "Random.H"
#include <tuple>
#include <chrono>

using namespace Foam;


template<typename T>
void testExclusivePrefixSum(T* result, T* input, label size, T zero_value)
{
    if (size <= 0) return;

    result[0] = zero_value;
    for (label i = 1; i < size; ++i)
    {
        result[i] = result[i - 1] + input[i - 1];
    }
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
        "Test for array prefix sum operation"
    );

    #include "setRootCase.H"
    #include "initDevice.H"
    #include "createMemoryPool.H"

    // include edge cases for warp,block and device wide prefix sum
    List<label> sizes = {0,1,8,32,33,64,65,800,1024,1025,4096,4097,6000,2*4096};

    using tuple_t =std::tuple<
        label,floatScalar,doubleScalar,
        Vector<label>,Vector<floatScalar>,Vector<doubleScalar>,
        Tensor<label>,Tensor<floatScalar>,Tensor<doubleScalar>
        >;
    tuple_t testType{};
    label Npassed = 0;
    const label numTypes = std::tuple_size<tuple_t>::value; // Number of types being tested
    const label Ntests = sizes.size()*numTypes;

    foamExecutor exec;

    // seed for rng
    auto now = std::chrono::high_resolution_clock::now();
    auto seed = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();

    auto testTypes = [&](auto typeTag)
    {
        using T = typename  decltype(typeTag)::type;
        Info << "Testing type: " << pTraits<T>::typeName << nl;

        for (label size : sizes)
        {
            Field<T> output(size);
            Field<T> input(size);
            Field<T> expected(size);

            // initialize input with random values

            Foam::Random rndGen(seed);
            for (size_t i = 0; i < size; ++i)
            {
                input[i] = randomValue<T>(rndGen);
            }

            // compute expected prefix sum
            testExclusivePrefixSum(expected.begin(), input.begin(), size, T(Zero));

            // compute output prefix sum
            exec.exclusivePrefixSum(output.begin(), input.begin(), size);

            // compare output with expected
            bool passed = true;
            const scalar componentTolerance = 1e-3;
            const scalar tolerance =
                componentTolerance * sqrt(pTraits<T>::nComponents);
            for (label i = 0; i < size; ++i)
            {
                if (mag(output[i] - expected[i]) > tolerance)
                {
                    passed = false;
                    Info << "\tMismatch at index " << i
                         << ": expected " << expected[i]
                         << ", got " << output[i]
                         << ", magDiff: " << mag(output[i] - expected[i])
                         << " , input: " << input[i] << nl;
                    break;
                }
            }

            Info << "\tprefix sum test, size: " << size << " on type: " <<
                pTraits<T>::typeName << " : " <<
                (passed ? "PASSED" : "FAILED") << nl;

            if (passed) ++Npassed;

        }

    };

    iterate_types(testTypes, std::move(testType));

    // Summary of all validation tests
    Info << "\n=== prefixSumTest execution summary: ===" << nl;
    Info << "\t testPassed: ["<< Npassed << " / " << Ntests << "]" << nl;

    return 0;
}

// ************************************************************************* //
