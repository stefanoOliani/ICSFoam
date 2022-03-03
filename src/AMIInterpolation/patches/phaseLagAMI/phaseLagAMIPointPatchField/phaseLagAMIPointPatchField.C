/*---------------------------------------------------------------------------*\

    ICSFoam: a library for Implicit Coupled Simulations in OpenFOAM
  
    Copyright (C) 2022  Stefano Oliani

    https://turbofe.it

-------------------------------------------------------------------------------
License
    This file is part of ICSFOAM.

    ICSFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ICSFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ICSFOAM.  If not, see <http://www.gnu.org/licenses/>.


Author
    Stefano Oliani
    Fluid Machinery Research Group, University of Ferrara, Italy
\*---------------------------------------------------------------------------*/

#include "phaseLagAMIPointPatchField.H"

#include "Swap.H"
#include "transformField.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::phaseLagAMIPointPatchField<Type>::phaseLagAMIPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(p, iF),
    phaseLagAMIPatch_(refCast<const phaseLagAMIPointPatch>(p)),
    ppiPtr_(nullptr),
    nbrPpiPtr_(nullptr)
{}


template<class Type>
Foam::phaseLagAMIPointPatchField<Type>::phaseLagAMIPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    coupledPointPatchField<Type>(p, iF, dict),
    phaseLagAMIPatch_(refCast<const phaseLagAMIPointPatch>(p, dict)),
    ppiPtr_(nullptr),
    nbrPpiPtr_(nullptr)
{
    if (!isType<phaseLagAMIPointPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "patch " << this->patch().index() << " not phaseLagAMI type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::phaseLagAMIPointPatchField<Type>::phaseLagAMIPointPatchField
(
    const phaseLagAMIPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    coupledPointPatchField<Type>(ptf, p, iF, mapper),
    phaseLagAMIPatch_(refCast<const phaseLagAMIPointPatch>(p)),
    ppiPtr_(nullptr),
    nbrPpiPtr_(nullptr)
{
    if (!isType<phaseLagAMIPointPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::phaseLagAMIPointPatchField<Type>::phaseLagAMIPointPatchField
(
    const phaseLagAMIPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(ptf, iF),
    phaseLagAMIPatch_(ptf.phaseLagAMIPatch_),
    ppiPtr_(nullptr),
    nbrPpiPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::phaseLagAMIPointPatchField<Type>::swapAddSeparated
(
    const Pstream::commsTypes,
    Field<Type>& pField
) const
{
    if (phaseLagAMIPatch_.phaseLagAMIPatch().owner())
    {
        // We inplace modify pField. To prevent the other side (which gets
        // evaluated at a later date) using already changed values we do
        // all swaps on the side that gets evaluated first.

        // Get neighbouring pointPatch
        const phaseLagAMIPointPatch& nbrPatch = phaseLagAMIPatch_.neighbPatch();

        // Get neighbouring pointPatchField
        const GeometricField<Type, pointPatchField, pointMesh>& fld =
            refCast<const GeometricField<Type, pointPatchField, pointMesh>>
            (
                this->internalField()
            );

        const phaseLagAMIPointPatchField<Type>& nbr =
            refCast<const phaseLagAMIPointPatchField<Type>>
            (
                fld.boundaryField()[nbrPatch.index()]
            );


        Field<Type> ptFld(this->patchInternalField(pField));
        Field<Type> nbrPtFld(nbr.patchInternalField(pField));


        if (doTransform())
        {
            const tensor& forwardT = this->forwardT()[0];
            const tensor& reverseT = this->reverseT()[0];

            transform(ptFld, reverseT, ptFld);
            transform(nbrPtFld, forwardT, nbrPtFld);
        }

        // convert point field to face field, AMI interpolate, then
        // face back to point
        {
            // add neighbour side contribution to owner
            Field<Type> nbrFcFld(nbrPpi().pointToFaceInterpolate(nbrPtFld));

            // interpolate to owner
            if (phaseLagAMIPatch_.phaseLagAMIPatch().applyLowWeightCorrection())
            {
                Field<Type> fcFld(ppi().pointToFaceInterpolate(ptFld));

                nbrFcFld =
                    phaseLagAMIPatch_.phaseLagAMIPatch().interpolate
                    (
                        nbrFcFld,
                        fcFld
                    );
            }
            else
            {
                nbrFcFld =
                    phaseLagAMIPatch_.phaseLagAMIPatch().interpolate(nbrFcFld);
            }

            // add to internal field
            this->addToInternalField
            (
                pField,
                ppi().faceToPointInterpolate(nbrFcFld)()
            );
        }

        {
            // add owner side contribution to neighbour
            Field<Type> fcFld(ppi().pointToFaceInterpolate(ptFld));

            // interpolate to neighbour
            if (phaseLagAMIPatch_.phaseLagAMIPatch().applyLowWeightCorrection())
            {
                Field<Type> nbrFcFld(nbrPpi().pointToFaceInterpolate(nbrPtFld));

                fcFld =
                    phaseLagAMIPatch_.phaseLagAMIPatch().neighbPatch().interpolate
                    (
                        fcFld,
                        nbrFcFld
                    );
            }
            else
            {
                fcFld =
                    phaseLagAMIPatch_.phaseLagAMIPatch().neighbPatch().interpolate
                    (
                        fcFld
                    );
            }

            // add to internal field
            nbr.addToInternalField
            (
                pField,
                nbrPpi().faceToPointInterpolate(fcFld)()
            );
        }
    }
}


// ************************************************************************* //
