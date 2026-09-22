#include "initParam.H"

internalField #codeStream
{
    code
    #{
        const IOdictionary& fieldDict =
            static_cast<const IOdictionary&>(dict);
        const fvMesh& mesh = refCast<const fvMesh>(fieldDict.db());
        scalarField field(mesh.nCells(), $p0);

        forAll(field, celli)
        {
            const vector& centre = mesh.C()[celli];
            const scalar x = centre.x();
            const scalar y = centre.y();
            const scalar z = centre.z();

            field[celli] =
                $p0
              + sqr($u0)/16
               *(cos(2*x/$L) + cos(2*y/$L))
               *(cos(2*z/$L) + 2);
        }

        field.writeEntry("", os);
    #};

    codeOptions
    #{
        -I$(LIB_SRC)/finiteVolume/lnInclude -I$(LIB_SRC)/meshTools/lnInclude
    #};

    codeInclude
    #{
        #include "fvCFD.H"
    #};
};
