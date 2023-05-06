/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

    Copyright (C) 2022 Stefano Oliani
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

Class
    Foam::phaseLagAMIPointPatchField

Description
    Cyclic AMI front and back plane patch field

SourceFiles
    phaseLagAMIPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef phaseLagAMIPointPatchField_H
#define phaseLagAMIPointPatchField_H

#include "phaseLagAMIPointPatch.H"
#include "coupledPointPatchField.H"
#include "PrimitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class phaseLagAMIPointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class phaseLagAMIPointPatchField
:
    public coupledPointPatchField<Type>
{
    // Private data

        //- Local reference cast into the phaseLagAMI patch
        const phaseLagAMIPointPatch& phaseLagAMIPatch_;

        //- Owner side patch interpolation pointer
        mutable autoPtr<PrimitivePatchInterpolation<primitivePatch>> ppiPtr_;

        //- Neighbour side patch interpolation pointer
        mutable autoPtr<PrimitivePatchInterpolation<primitivePatch>>
            nbrPpiPtr_;


    // Private Member Functions

        //- Owner side patch interpolation
        const PrimitivePatchInterpolation<primitivePatch>& ppi() const
        {
            if (!ppiPtr_.valid())
            {
                ppiPtr_.reset
                (
                    new PrimitivePatchInterpolation<primitivePatch>
                    (
                        phaseLagAMIPatch_.phaseLagAMIPatch()
                    )
                );
            }

            return *ppiPtr_;
        }

        //- Neighbour side patch interpolation
        const PrimitivePatchInterpolation<primitivePatch>& nbrPpi() const
        {
            if (!nbrPpiPtr_.valid())
            {
                nbrPpiPtr_.reset
                (
                    new PrimitivePatchInterpolation<primitivePatch>
                    (
                        phaseLagAMIPatch_.phaseLagAMIPatch().neighbPatch()
                    )
                );
            }

            return *nbrPpiPtr_;
        }


public:

    //- Runtime type information
    TypeName(phaseLagAMIPointPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        phaseLagAMIPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        phaseLagAMIPointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        phaseLagAMIPointPatchField
        (
            const phaseLagAMIPointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<Type>> clone() const
        {
            return autoPtr<pointPatchField<Type>>
            (
                new phaseLagAMIPointPatchField<Type>
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        phaseLagAMIPointPatchField
        (
            const phaseLagAMIPointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<Type>>
            (
                new phaseLagAMIPointPatchField<Type>
                (
                    *this, iF
                )
            );
        }


    // Member functions

        // Constraint handling

            //- Return the constraint type this pointPatchField implements
            virtual const word& constraintType() const
            {
                return phaseLagAMIPointPatch::typeName;
            }


        // Cyclic AMI coupled interface functions

            //- Does the patch field perform the transformation
            virtual bool doTransform() const
            {
                return
                   !(
                        phaseLagAMIPatch_.parallel()
                     || pTraits<Type>::rank == 0
                    );
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return phaseLagAMIPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return phaseLagAMIPatch_.reverseT();
            }


        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            )
            {}

            //- Complete swap of patch point values and add to local values
            virtual void swapAddSeparated
            (
                const Pstream::commsTypes commsType,
                Field<Type>&
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "phaseLagAMIPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
