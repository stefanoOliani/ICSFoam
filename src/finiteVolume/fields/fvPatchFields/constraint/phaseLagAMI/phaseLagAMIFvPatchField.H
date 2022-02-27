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

#ifndef phaseLagAMIFvPatchField_H
#define phaseLagAMIFvPatchField_H

#include "coupledFvPatchField.H"
#include "phaseLagAMILduInterfaceField.H"
#include "phaseLagAMIFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class phaseLagAMIFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class phaseLagAMIFvPatchField
:
    virtual public phaseLagAMILduInterfaceField,
    public coupledFvPatchField<Type>
{
    // Private data

        //- Local reference cast into the cyclic patch
        const phaseLagAMIFvPatch& phaseLagAMIPatch_;

    // Private Member Functions

        //- Return neighbour side field given internal fields
        template<class Type2>
        tmp<Field<Type2>> neighbourSideField
        (
            const Field<Type2>&
        ) const;


public:

    //- Runtime type information
    TypeName(phaseLagAMIFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        phaseLagAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        phaseLagAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given phaseLagAMIFvPatchField onto a new patch
        phaseLagAMIFvPatchField
        (
            const phaseLagAMIFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        phaseLagAMIFvPatchField(const phaseLagAMIFvPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new phaseLagAMIFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        phaseLagAMIFvPatchField
        (
            const phaseLagAMIFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new phaseLagAMIFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return local reference cast into the cyclic AMI patch
            const phaseLagAMIFvPatch& phaseLagAMIPatch() const
            {
                return phaseLagAMIPatch_;
            }


        // Evaluation functions

            //- Return true if coupled. Note that the underlying patch
            //  is not coupled() - the points don't align.
            virtual bool coupled() const;

            //- Return neighbour coupled internal cell data
            virtual tmp<Field<Type>> patchNeighbourField() const;

            //- Return reference to neighbour patchField
            const phaseLagAMIFvPatchField<Type>& neighbourPatchField() const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                solveScalarField& result,
                const bool add,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>&,
                const bool add,
                const Field<Type>&,
                const scalarField&,
                const Pstream::commsTypes commsType
            ) const;


        // Cyclic AMI coupled interface functions

            //- Does the patch field perform the transformation
            virtual bool doTransform() const
            {
                return
                    !(phaseLagAMIPatch_.parallel() || pTraits<Type>::rank == 0);
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

            //- Return rank of component for transform
            virtual int rank() const
            {
                return pTraits<Type>::rank;
            }


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


//- Update result field based on interface functionality
template<>
void phaseLagAMIFvPatchField<scalar>::updateInterfaceMatrix
(
	solveScalarField& result,
	const bool add,
	const solveScalarField& psiInternal,
	const scalarField& coeffs,
	const direction cmpt,
	const Pstream::commsTypes commsType
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "phaseLagAMIFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
