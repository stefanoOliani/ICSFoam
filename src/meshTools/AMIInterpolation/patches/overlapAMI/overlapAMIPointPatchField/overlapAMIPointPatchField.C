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

#include "overlapAMIPointPatchField.H"

#include "Swap.H"
#include "transformField.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::overlapAMIPointPatchField<Type>::overlapAMIPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(p, iF),
    overlapAMIPatch_(refCast<const overlapAMIPointPatch>(p)),
    ppiPtr_(nullptr),
    nbrPpiPtr_(nullptr)
{}


template<class Type>
Foam::overlapAMIPointPatchField<Type>::overlapAMIPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    coupledPointPatchField<Type>(p, iF, dict),
    overlapAMIPatch_(refCast<const overlapAMIPointPatch>(p, dict)),
    ppiPtr_(nullptr),
    nbrPpiPtr_(nullptr)
{
    if (!isType<overlapAMIPointPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "patch " << this->patch().index() << " not overlapAMI type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::overlapAMIPointPatchField<Type>::overlapAMIPointPatchField
(
    const overlapAMIPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    coupledPointPatchField<Type>(ptf, p, iF, mapper),
    overlapAMIPatch_(refCast<const overlapAMIPointPatch>(p)),
    ppiPtr_(nullptr),
    nbrPpiPtr_(nullptr)
{
    if (!isType<overlapAMIPointPatch>(this->patch()))
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
Foam::overlapAMIPointPatchField<Type>::overlapAMIPointPatchField
(
    const overlapAMIPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    coupledPointPatchField<Type>(ptf, iF),
    overlapAMIPatch_(ptf.overlapAMIPatch_),
    ppiPtr_(nullptr),
    nbrPpiPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::overlapAMIPointPatchField<Type>::swapAddSeparated
(
    const Pstream::commsTypes,
    Field<Type>& pField
) const
{
    if (overlapAMIPatch_.overlapAMIPatch().owner())
    {
        // We inplace modify pField. To prevent the other side (which gets
        // evaluated at a later date) using already changed values we do
        // all swaps on the side that gets evaluated first.

        // Get neighbouring pointPatch
        const overlapAMIPointPatch& nbrPatch = overlapAMIPatch_.neighbPatch();

        // Get neighbouring pointPatchField
        const GeometricField<Type, pointPatchField, pointMesh>& fld =
            refCast<const GeometricField<Type, pointPatchField, pointMesh>>
            (
                this->internalField()
            );

        const overlapAMIPointPatchField<Type>& nbr =
            refCast<const overlapAMIPointPatchField<Type>>
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
            if (overlapAMIPatch_.overlapAMIPatch().applyLowWeightCorrection())
            {
                Field<Type> fcFld(ppi().pointToFaceInterpolate(ptFld));

                nbrFcFld =
                    overlapAMIPatch_.overlapAMIPatch().interpolate
                    (
                        nbrFcFld,
                        fcFld
                    );
            }
            else
            {
                nbrFcFld =
                    overlapAMIPatch_.overlapAMIPatch().interpolate(nbrFcFld);
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
            if (overlapAMIPatch_.overlapAMIPatch().applyLowWeightCorrection())
            {
                Field<Type> nbrFcFld(nbrPpi().pointToFaceInterpolate(nbrPtFld));

                fcFld =
                    overlapAMIPatch_.overlapAMIPatch().neighbPatch().interpolate
                    (
                        fcFld,
                        nbrFcFld
                    );
            }
            else
            {
                fcFld =
                    overlapAMIPatch_.overlapAMIPatch().neighbPatch().interpolate
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
