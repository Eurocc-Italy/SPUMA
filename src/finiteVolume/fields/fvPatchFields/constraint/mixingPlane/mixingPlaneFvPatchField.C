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
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMatrix.H"
#include "mixingPlaneFvPatchField.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "fvsPatchField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::mixingPlaneFvPatchField<Type>::coupled() const
{
    // return mixingPlanePatch_.coupled();
    return true;
}

template<class Type>
void mixingPlaneFvPatchField<Type>::readMixingType() const
{
    dictionary schemesDict;

    if(this->patch().boundaryMesh().mesh().schemesDict().found("mixingPlane"))
    {
        schemesDict = this->patch().boundaryMesh().mesh().schemesDict().subDict
        (
            "mixingPlane"
        );
    }
    else
    {
        // Probably running from reconstructPar, so an empty dictionary
        // Reading from the case's system directory
        IOdictionary dict
        (
            IOobject
            (
                "fvSchemes",
                this->patch().boundaryMesh().mesh().time().caseSystem(),
                this->patch().boundaryMesh().mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        );

        schemesDict = dict.subDict
        (
            "mixingPlane"
        );
    }

    // Try reading field type
    word fieldName = this->internalField().name();

    if (schemesDict.found(fieldName))
    {
        mixing_ = mixingPlaneInterpolation::mixingTypeNames_.read
        (
            schemesDict.lookup(fieldName)
        );
    }
    else if (schemesDict.found("default"))
    {
        mixing_ = mixingPlaneInterpolation::mixingTypeNames_.read
        (
            schemesDict.lookup("default")
        );
    }
    else
    {
        FatalIOErrorIn
        (
            "void mixingPlaneFvPatchField<Type>::readMixingType() const",
            schemesDict
        )   << "Cannot find mixing type for field "
            <<  this->internalField().name() << nl
            << "Please specify in fvSchemes in mixingPlane, "
            << "under field name " << nl
            << "Available types are "
            << mixingPlaneInterpolation::mixingTypeNames_.toc()
            << abort(FatalIOError);
    }
}


template<class Type>
void mixingPlaneFvPatchField<Type>::calcFluxMask() const
{
    // Find the flux field and calculate flux mask
    if (mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING)
    {
        if (!this->db().objectRegistry::found(phiName_))
        {
            InfoIn
            (
                "void mixingPlaneFvPatchField<Type>::calcFluxMask()"
                ""
            )   << "Flux not found for flux averaging mixing plane on "
                << this->patch().name() << " for field "
                << this->internalField().name()
                << ".  Flux field: " << phiName_
                << endl;

            fluxMask_.setSize(mixingPlanePatch_.nProfileBands(), 0);
        }
        else
        {
            // Get flux field.  Used in decision for which faces
            // see the average and for flux-based weights
            const scalarField& phip = this->patch().template lookupPatchField<surfaceScalarField>
            (
                phiName_
            );

            if (mixingPlanePatch_.master())
            {
                // Calculate flux masks.  Master fluxes are used for the mask
                // Flux mask behaves like valueFraction:
                //     1 = use interpolated value
                //     0 = use patch internal field

                // Note: only master calculates the flux mask, as it is
                // contained on the profile.  Shadow will use the flux mask
                // from master.  HJ, 14/Feb/2013
                fluxMask_ = neg(mixingPlanePatch_.toProfile(phip));
            }
            else
            {
                // On the shadow side, flux mask is opposite from master
                fluxMask_ = 1 - shadowPatchField().fluxMask_;
            }

            // Calculate flux weights via circumferential average
            // Note: for zero flux, weights must be uniform and not zero
            // HJ, 4/Feb/2013

            // Note: flux weights exist on both master and slave
            scalar maxPhiP = 0;

            if (!phip.empty())
            {
                maxPhiP = max(phip);
            }

            reduce(maxPhiP, maxOp<scalar>());

            if (maxPhiP > SMALL)
            {
                fluxWeights_ = Foam::max(phip, scalar(0));
                fluxWeights_ /=
                    mixingPlanePatch_.circumferentialAverage(fluxWeights_)
                  + SMALL;
           }
            else
            {
                fluxWeights_ = 1;
            }
        }
    }
}


template<class Type>
const scalarField& mixingPlaneFvPatchField<Type>::fluxMask() const
{
    if (!this->updated())
    {
        // Force recalculation of flux masks
        calcFluxMask();
    }

    return fluxMask_;
}


template<class Type>
const scalarField& mixingPlaneFvPatchField<Type>::fluxWeights() const
{
    if (!this->updated())
    {
        // Force recalculation of flux masks
        calcFluxMask();
    }

    return fluxWeights_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(p)),
    mixing_(mixingPlaneInterpolation::MIXING_UNKNOWN),
    phiName_("phi"),
    fluxMask_(),
    fluxWeights_(p.size(), 0)
{
    // Cannot read mixing type here: internal field may not be set
    // HJ, 3/Jun/2015
}


template<class Type>
mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixingPlaneLduInterfaceField(),
    coupledFvPatchField<Type>(p, iF, dict, false),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(p)),
    mixing_(mixingPlaneInterpolation::MIXING_UNKNOWN),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    fluxMask_(),
    fluxWeights_(p.size(), 0)
{
    if (!isType<mixingPlaneFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not mixingPlane type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }

    // Read mixing type
    this->readMixingType();

    if (dict.found("value"))
    {
#if OPENFOAM < 2302
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
#else
        this->readValueEntry(dict);
#endif
    }
    else
    {
        // Grab the internal value for initialisation.
        fvPatchField<Type>::operator=(this->patchInternalField()());
    }
}


template<class Type>
mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const mixingPlaneFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(p)),
    mixing_(ptf.mixing_),
    phiName_(ptf.phiName_),
    fluxMask_(),
    fluxWeights_(ptf.fluxWeights_, mapper)
{
    if (!isType<mixingPlaneFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField\n"
            "(\n"
            "    const mixingPlaneFvPatchField<Type>& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const mixingPlaneFvPatchField<Type>& ptf
)
:
    mixingPlaneLduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(ptf.patch())),
    mixing_(ptf.mixing_),
    phiName_(ptf.phiName_),
    fluxMask_(),
    fluxWeights_(ptf.fluxWeights_)
{}


template<class Type>
mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const mixingPlaneFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixingPlaneLduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(ptf.patch())),
    mixing_(ptf.mixing_),
    phiName_(ptf.phiName_),
    fluxMask_(),
    fluxWeights_(ptf.fluxWeights_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void mixingPlaneFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fluxMask_.clear();
    fluxWeights_.autoMap(m);
}


template<class Type>
void mixingPlaneFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    fluxMask_.clear();

    const mixingPlaneFvPatchField<Type>& tiptf =
        refCast<const mixingPlaneFvPatchField<Type> >(ptf);

    fluxWeights_.rmap(tiptf.fluxWeights_, addr);
}


// Return shadow field
template<class Type>
const mixingPlaneFvPatchField<Type>&
mixingPlaneFvPatchField<Type>::shadowPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->internalField()
        );

    return refCast<const mixingPlaneFvPatchField<Type> >
    (
        fld.boundaryField()[mixingPlanePatch_.shadowIndex()]
    );
}


// Return shadow field
template<class Type>
const mixingPlaneInterpolation::mixingType&
mixingPlaneFvPatchField<Type>::mixing() const
{
    // If mixing type is unknown, read it
    if (this->mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
    {
        this->readMixingType();
    }

    return mixing_;
}


template<class Type>
void mixingPlaneFvPatchField<Type>::updateCoeffs()
{
    // Read mixing type
    readMixingType();

    // Force recalculation of flux masks
    calcFluxMask();

    fvPatchField<Type>::updateCoeffs();
}


template<class Type>
tmp<Field<Type>> mixingPlaneFvPatchField<Type>::patchNeighbourField() const
{
    // Get shadow patch internalField field
    const Field<Type> sField(this->shadowPatchField().patchInternalField());

    // If mixing type is unknown, read it
    if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
    {
        this->readMixingType();
    }

    if (mixing_ == mixingPlaneInterpolation::AREA_AVERAGING)
    {
        // Area-weighted averaging
        return mixingPlanePatch_.interpolate(sField);
    }
    else if (mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING)
    {
        // Flux averaging
        // - for outgoing flux, use zero gradient condition
        // - for incoming flux, use interpolated flux-weighted value

        const scalarField& mask = fluxMask();

        const scalarField& shadowFluxWeights =
            shadowPatchField().fluxWeights();

        // For outgoing flux, the value is identical to internal value
        // For incoming flux, calculate the average value of the
        // flux-weight shadow values coming out
        return
            mixingPlanePatch_.fromProfile
            (
                mask*mixingPlanePatch_.shadow().toProfile
                (
                    sField*shadowFluxWeights
                )
            )
          + mixingPlanePatch_.fromProfile(1 - mask)*this->patchInternalField();
    }
    else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
    {
        return this->patchInternalField();
    }
    else
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > mixingPlaneFvPatchField<Type>::"
            "patchNeighbourField() const"
        )   << "Unknown mixing type for patch " << this->patch().name()
            << " for field "
            <<  this->internalField().name()
            << abort(FatalError);
    }

    // Dummy return to keep compiler happy
    return this->patchInternalField();
}

template<class Type>
void Foam::mixingPlaneFvPatchField<Type>::patchNeighbourField
(
    UList<Type>& pnf
) const
{
    // checkCommunicator = false
    auto tpnf = this->patchNeighbourField();
    pnf.deepCopy(tpnf());
}


template<class Type>
tmp<scalarField>
mixingPlaneFvPatchField<Type>::untransformedInterpolate
(
    const direction cmpt
) const
{
    const Field<Type>& iField = this->internalField();

    // Get shadow face-cells and assemble shadow field
    const labelUList& sfc = mixingPlanePatch_.shadow().faceCells();

    scalarField sField(sfc.size());

    forAll (sField, i)
    {
        sField[i] = component(iField[sfc[i]], cmpt);
    }

    // if (mixing_ == mixingPlaneInterpolation::AREA_AVERAGING)
    {
        // Area-weighted averaging
        return mixingPlanePatch_.interpolate(sField);
    }
    // else if (mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING)
    // {
    //     // Flux averaging
    //     // - for outgoing flux, use zero gradient condition
    //     // - for incoming flux, use interpolated flux-weighted value

    //     //HJ, HERE:!!!
    // }
    // else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
    // {
    //     return this->patchInternalField()().component(cmpt);
    // }
    // else
    // {
    //     FatalErrorIn
    //     (
    //         "tmp<scalarField> mixingPlaneFvPatchField<Type>::"
    //         "untransformedInterpolate(const direction cmpt) const"
    //     )   << "Unknown mixing type for patch " << this->patch().name()
    //         << " for field "
    //         <<  this->internalField().name()
    //         << abort(FatalError);
    // }
}


template<class Type>
void mixingPlaneFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Evaluation is identical with or without flux averaging
    // as the masking and weighting is taken into account in
    // patchNeighbourField.  HJ, 4/Feb/2013

    const scalarField& w = this->patch().weights();

    // If mixing type is unknown, read it
    if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
    {
        this->readMixingType();
    }

    if
    (
        mixing_ == mixingPlaneInterpolation::AREA_AVERAGING
     || mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING
    )
    {
        Field<Type>::operator=
        (
            w*this->mixingPlanePatch_.circumferentialAverage
            (
                this->patchInternalField()
            )
          + (1 - w)*this->patchNeighbourField()
        );
    }
    else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
    {
        Field<Type>::operator=(this->patchInternalField());
    }
    else
    {
        FatalErrorIn
        (
            "void mixingPlaneFvPatchField<Type>::initEvaluate\n"
            "(\n"
            "    const Pstream::commsTypes commsType\n"
            ")"
        )   << "Unknown mixing type for patch " << this->patch().name()
            << " for field "
            <<  this->internalField().name()
            << abort(FatalError);
    }
}


template<class Type>
void mixingPlaneFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsTypes
)
{
    coupledFvPatchField<Type>::evaluate(commsTypes);
}


template<class Type>
tmp<Field<Type> > mixingPlaneFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    // If mixing type is unknown, read it
    if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
    {
        this->readMixingType();
    }

    if (mixing_ == mixingPlaneInterpolation::AREA_AVERAGING)
    {
        return pTraits<Type>::one*w;
    }
    else if (mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING)
    {
        // Flux averaging.
        // When flux averaging indicates outgoing flux,
        // use zero gradient condition: evaluating as internal field
        // When flux averaging indicates incoming flux,
        // use normal interpolation
        // Note that each face may have elements in multiple bands and
        // therefore interpolation from profile to patch is required
        // HJ, 11/Feb/2013
        const scalarField& mask = fluxMask();

        scalarField oneMFluxMask(mixingPlanePatch_.fromProfile(1 - mask));

        return pTraits<Type>::one*oneMFluxMask;
    }
    else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
    {
        return tmp<Field<Type> >
        (
            new Field<Type>(this->size(), pTraits<Type>::one)
        );
    }
    else
    {
        FatalErrorIn
        (
            "void mixingPlaneFvPatchField<Type>::valueInternalCoeffs\n"
            "(\n"
            "    const tmp<scalarField>& w\n"
            ") const"
        )   << "Unknown mixing type for patch " << this->patch().name()
            << " for field "
            <<  this->internalField().name()
            << abort(FatalError);
    }

    // Dummy return to keep compiler happy
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}


template<class Type>
tmp<Field<Type> > mixingPlaneFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    // If mixing type is unknown, read it
    if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
    {
        this->readMixingType();
    }

    if (mixing_ == mixingPlaneInterpolation::AREA_AVERAGING)
    {
        return pTraits<Type>::one*(1.0 - w);
    }
    else if (mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING)
    {
        // Flux averaging
        const scalarField& mask = fluxMask();

        const scalarField fluxMask(mixingPlanePatch_.fromProfile(mask));

        return pTraits<Type>::one*fluxMask;
    }
    else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
    {
        return  tmp<Field<Type> >
        (
            new Field<Type>(this->size(), pTraits<Type>::zero)
        );
    }
    else
    {
        FatalErrorIn
        (
            "void mixingPlaneFvPatchField<Type>::valueBoundaryCoeffs\n"
            "(\n"
            "    const tmp<scalarField>& w\n"
            ") const"
        )   << "Unknown mixing type for patch " << this->patch().name()
            << " for field "
            <<  this->internalField().name()
            << abort(FatalError);
    }

    // Dummy return to keep compiler happy
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}


// template<class Type>
// tmp<Field<Type> >
// mixingPlaneFvPatchField<Type>::gradientInternalCoeffs() const
// {
//     // If mixing type is unknown, read it
//     if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
//     {
//         this->readMixingType();
//     }
//
//     if
//     (
//         mixing_ == mixingPlaneInterpolation::AREA_AVERAGING
//      || mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING
//     )
//     {
//         return -pTraits<Type>::one*this->patch().deltaCoeffs();
//     }
//     else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
//     {
//         return tmp<Field<Type> >
//         (
//             new Field<Type>(this->size(), pTraits<Type>::zero)
//         );
//     }
//     else
//     {
//         FatalErrorIn
//         (
//             "void mixingPlaneFvPatchField<Type>::gradientInternalCoeffs"
//             "gradientInternalCoeffs() const"
//         )   << "Unknown mixing type for patch " << this->patch().name()
//             << " for field "
//             <<  this->internalField().name()
//             << abort(FatalError);
//     }

//     // Dummy return to keep compiler happy
//     return tmp<Field<Type> >
//     (
//         new Field<Type>(this->size(), pTraits<Type>::zero)
//     );
// }


// template<class Type>
// tmp<Field<Type> >
// mixingPlaneFvPatchField<Type>::gradientBoundaryCoeffs() const
// {
//     // If mixing type is unknown, read it
//     if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
//     {
//         this->readMixingType();
//     }
//
//     if
//     (
//         mixing_ == mixingPlaneInterpolation::AREA_AVERAGING
//      || mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING
//     )
//     {
//         -this->gradientInternalCoeffs();
//     }
//     else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
//     {
//         return tmp<Field<Type> >
//         (
//             new Field<Type>(this->size(), pTraits<Type>::zero)
//         );
//     }
//     else
//     {
//         FatalErrorIn
//         (
//             "void mixingPlaneFvPatchField<Type>::gradientBoundaryCoeffs"
//             "gradientBoundaryCoeffs() const"
//         )   << "Unknown mixing type for patch " << this->patch().name()
//             << " for field "
//             <<  this->internalField().name()
//             << abort(FatalError);
//     }

//     // Dummy return to keep compiler happy
//     return tmp<Field<Type> >
//     (
//         new Field<Type>(this->size(), pTraits<Type>::zero)
//     );
// }


template<class Type>
void Foam::mixingPlaneFvPatchField<Type>::initInterfaceMatrixUpdate
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label interfacei,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
}

template<class Type>
void Foam::mixingPlaneFvPatchField<Type>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    // Communication is allowed either before or after processor
    // patch comms.  HJ, 11/Jul/2011

    // If mixing type is unknown, read it
    if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
    {
        this->readMixingType();
    }

    if
    (
        mixing_ == mixingPlaneInterpolation::AREA_AVERAGING
     || mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING
    )
    {
        // Get shadow face-cells and assemble shadow field
        const labelUList& sfc = mixingPlanePatch_.shadow().faceCells();

        scalarField sField(sfc.size());

        forAll (sField, i)
        {
            sField[i] = psiInternal[sfc[i]];
        }

        // Get local faceCells
        const labelUList& fc = mixingPlanePatch_.faceCells();

        scalarField pnf(mixingPlanePatch_.interpolate(sField));

        // Multiply the field by coefficients and add into the result
        this->addToInternalField(result, !add, fc, coeffs, pnf);
    }
    else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
    {
        // Do nothing
    }
    else
    {
        FatalErrorIn
        (
            "void mixingPlaneFvPatchField<Type>::updateInterfaceMatrixUpdate\n"
            "(\n"
            "    const scalarField& psiInternal,\n"
            "    scalarField& result,\n"
            "    const lduMatrix&,\n"
            "    const scalarField& coeffs,\n"
            "    const direction cmpt,\n"
            "    const Pstream::commsTypes commsType,\n"
            "    const bool switchToLhs\n"
            ") const"
        )   << "Unknown mixing type for patch " << this->patch().name()
            << " for field "
            <<  this->internalField().name()
            << abort(FatalError);
    }
}

template<class Type>
void Foam::mixingPlaneFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{
    Info << "DEBUG13 need to check: " << this->patch().name() << " : " << this->internalField().name() << endl;
    // Communication is allowed either before or after processor
    // patch comms.  HJ, 11/Jul/2011

    // If mixing type is unknown, read it
    if (mixing_ == mixingPlaneInterpolation::MIXING_UNKNOWN)
    {
        this->readMixingType();
    }

    if
    (
        mixing_ == mixingPlaneInterpolation::AREA_AVERAGING
     || mixing_ == mixingPlaneInterpolation::FLUX_AVERAGING
    )
    {
        // Get shadow face-cells and assemble shadow field
        const labelUList& sfc = mixingPlanePatch_.shadow().faceCells();

        Field<Type> sField(sfc.size());

        forAll (sField, i)
        {
            sField[i] = psiInternal[sfc[i]];
        }

        // Get local faceCells
        const labelUList& fc = mixingPlanePatch_.faceCells();
        const Field<Type> pnf(mixingPlanePatch_.interpolate(sField));

        // Multiply the field by coefficients and add into the result
        this->addToInternalField(result, !add, fc, coeffs, pnf);

    }
    else if (mixing_ == mixingPlaneInterpolation::ZERO_GRADIENT)
    {
        // Do nothing
    }
    else
    {
        FatalErrorIn
        (
            "void mixingPlaneFvPatchField<Type>::updateInterfaceMatrixUpdate\n"
            "(\n"
            "    const scalarField& psiInternal,\n"
            "    scalarField& result,\n"
            "    const lduMatrix&,\n"
            "    const scalarField& coeffs,\n"
            "    const direction cmpt,\n"
            "    const Pstream::commsTypes commsType,\n"
            "    const bool switchToLhs\n"
            ") const"
        )   << "Unknown mixing type for patch " << this->patch().name()
            << " for field "
            <<  this->internalField().name()
            << abort(FatalError);
    }
}

template<class Type>
void Foam::mixingPlaneFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix,
    const label mat,
    const direction cmpt
)
{
    Info << "DEBUG14: need to check: " << this->patch().name() << " : " << this->internalField().name() << endl;

    if (this->mixingPlanePatch_.owner())
    {
        label index = this->patch().index();

        const label globalPatchID =
            matrix.lduMeshAssembly().patchLocalToGlobalMap()[mat][index];

        const Field<scalar> intCoeffsCmpt
        (
            matrix.internalCoeffs()[globalPatchID].component(cmpt)
        );

        const Field<scalar> boundCoeffsCmpt
        (
            matrix.boundaryCoeffs()[globalPatchID].component(cmpt)
        );

        tmp<Field<scalar>> tintCoeffs(coeffs(matrix, intCoeffsCmpt, mat));
        tmp<Field<scalar>> tbndCoeffs(coeffs(matrix, boundCoeffsCmpt, mat));
        const Field<scalar>& intCoeffs = tintCoeffs.ref();
        const Field<scalar>& bndCoeffs = tbndCoeffs.ref();

        const labelUList& u = matrix.lduAddr().upperAddr();
        const labelUList& l = matrix.lduAddr().lowerAddr();

        const labelList& faceMap =
            matrix.lduMeshAssembly().faceBoundMap()[mat][index];

        forAll (faceMap, faceI)
        {
            label globalFaceI = faceMap[faceI];

            const scalar boundCorr = -bndCoeffs[faceI];
            const scalar intCorr = -intCoeffs[faceI];

#ifdef _FOLLOW_CYCLICFVPATCHFIELD_COMPUTATION
            // This is taken from cyclicFvPatchField.C
            matrix.upper()[globalFaceI] += boundCorr;
            matrix.diag()[u[globalFaceI]] -= boundCorr; // <- those 2 lines?
            matrix.diag()[l[globalFaceI]] -= intCorr;   // <-

#else
            // This is taken from cyclicACMIFvPatchField.C
            matrix.upper()[globalFaceI] += boundCorr;
            matrix.diag()[u[globalFaceI]] -= intCorr;   // <- or those 2??
            matrix.diag()[l[globalFaceI]] -= boundCorr; // <-

#endif
            if (matrix.asymmetric())
            {
                matrix.lower()[globalFaceI] += intCorr;
            }
        }

        if (matrix.psi(mat).mesh().fluxRequired(this->internalField().name()))
        {
            matrix.internalCoeffs().set
            (
                globalPatchID, intCoeffs*pTraits<Type>::one
            );
            matrix.boundaryCoeffs().set
            (
                globalPatchID, bndCoeffs*pTraits<Type>::one
            );

            const label nbrPathID = mixingPlanePatch_.shadowIndex();

            const label nbrGlobalPatchID =
                matrix.lduMeshAssembly().patchLocalToGlobalMap()[mat][nbrPathID];

            matrix.internalCoeffs().set
            (
                nbrGlobalPatchID, intCoeffs*pTraits<Type>::one

            );
            matrix.boundaryCoeffs().set
            (
                nbrGlobalPatchID, bndCoeffs*pTraits<Type>::one
            );
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Foam::scalar>>
Foam::mixingPlaneFvPatchField<Type>::coeffs
(
    fvMatrix<Type>& matrix,
    const Field<scalar>& coeffs,
    const label mat
) const
{
    Info << "DEBUG15: need to check: " << this->patch().name() << " : " << this->internalField().name() << endl;
    const label index(this->patch().index());

    const label nSubFaces
    (
        matrix.lduMeshAssembly().cellBoundMap()[mat][index].size()
    );

    Field<scalar> mapCoeffs(nSubFaces, Zero);

    const scalarListList& srcWeight =
        mixingPlanePatch_.mixingPlanePatch().patchToPatch().slaveProfileToPatchWeights();

    label subFaceI = 0;
    forAll(*this, faceI)
    {
        const scalarList& w = srcWeight[faceI];
        for(label i=0; i<w.size(); i++)
        {
            const label localFaceId =
                matrix.lduMeshAssembly().facePatchFaceMap()[mat][index][subFaceI];
            mapCoeffs[subFaceI] = w[i]*coeffs[localFaceId];
            subFaceI++;
        }
    }

    return tmp<Field<scalar>>(new Field<scalar>(mapCoeffs));
}

template<class Type>
void mixingPlaneFvPatchField<Type>::writeAveragedField(Ostream& os) const
{
  // Create a temporary clone from this patchField
  tmp<mixingPlaneFvPatchField<Type> > tmpF(new mixingPlaneFvPatchField<Type>(*this));
  mixingPlaneFvPatchField<Type>& avgF = tmpF.ref();

  // Re-evaluate the mixingPlane field for this copy based on the selected averaging scheme
  // (areaAveraging or fluxAveraging)
  avgF.initEvaluate(UPstream::defaultCommsType);

  // Write to file
  avgF.writeEntry("value", os);
}

template<class Type>
void mixingPlaneFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    os.writeEntryIfDifferent("phi", word("phi"), phiName_);

    // Recompute field's average before writing to disk?
    if
    (
        mixingPlanePatch_.writeAveragedPatchField() &&
        this->db().objectRegistry::found(phiName_)      // true if solver is active
    )
    {
        this->writeAveragedField(os);
    }
    else
    {
        this->writeEntry("value", os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
