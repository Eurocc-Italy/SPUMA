#include "tensor4.H"
#include "vector4.H"
#include "IOstreams.H"

using namespace Foam;

int main()
{
    const tensor4 matrix
    (
        4, 1, 0, 0,
        0, 3, 1, 0,
        0, 0, 2, 1,
        1, 0, 0, 5
    );
    const vector4 expected(1, 2, 3, 4);
    const vector4 rhs(6, 9, 10, 21);
    const vector4 product = dot(matrix, expected);
    const vector4 operatorProduct = matrix & expected;
    const vector4 recovered = solve(matrix, rhs);

    for (direction component = 0; component < 4; ++component)
    {
        if
        (
            mag(product[component] - rhs[component]) > SMALL
         || mag(operatorProduct[component] - rhs[component]) > SMALL
         || mag(recovered[component] - expected[component]) > SMALL
        )
        {
            FatalErrorInFunction
                << "Tensor4 primitive algebra test failed" << exit(FatalError);
        }
    }

    Info<< "Tensor4 primitive algebra test passed" << endl;
    return 0;
}
