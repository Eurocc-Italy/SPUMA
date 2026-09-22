#include "initParam.H"

internalField #codeStream
{
    code
    #{
        const IOdictionary& fieldDict =
            static_cast<const IOdictionary&>(dict);
        const fvMesh& mesh = refCast<const fvMesh>(fieldDict.db());
        vectorField field(mesh.nCells(), vector::zero);

        forAll(field, celli)
        {
            const vector& centre = mesh.C()[celli];
            const scalar x = centre.x();
            const scalar y = centre.y();
            const scalar z = centre.z();

            field[celli] = vector
            (
                 $u0*sin(x/$L)*cos(y/$L)*cos(z/$L),
                -$u0*cos(x/$L)*sin(y/$L)*cos(z/$L),
                 0
            );
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
