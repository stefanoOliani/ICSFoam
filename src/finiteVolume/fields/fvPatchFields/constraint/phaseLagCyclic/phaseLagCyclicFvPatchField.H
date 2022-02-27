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

#ifndef phaseLagCyclicFvPatchField_H
#define phaseLagCyclicFvPatchField_H

#include "cyclicFvPatchField.H"
#include "cyclicLduInterfaceField.H"
#include "phaseLagCyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class phaseLagCyclicFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class phaseLagCyclicFvPatchField
:
	public cyclicFvPatchField<Type>
{

    // Private data

		//- Sub time level number
		label subTimeLevel_;

		//- Local reference cast into the cyclic patch
		const phaseLagCyclicFvPatch& phaseLagCyclicPatch_;


public:

    //- Runtime type information
    TypeName("phaseLagCyclic");

    // Constructors

        //- Construct from patch and internal field
        phaseLagCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        phaseLagCyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given phaseLagCyclicFvPatchField onto a
        //  new patch
        phaseLagCyclicFvPatchField
        (
            const phaseLagCyclicFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        phaseLagCyclicFvPatchField
        (
            const phaseLagCyclicFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new phaseLagCyclicFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        phaseLagCyclicFvPatchField
        (
            const phaseLagCyclicFvPatchField<Type>&,
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
                new phaseLagCyclicFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

			//- Return the interface type
			virtual const word& interfaceFieldType() const
			{
				return cyclicFvPatchField<Type>::type();
			}

            //- Return local reference cast into the cyclic patch
            const phaseLagCyclicFvPatch& phaseLagCyclicPatch() const
            {
                return phaseLagCyclicPatch_;
            }

            //- Return the "jump" across the patch
            virtual tmp<Field<Type>> phaseLaggedField() const;

		// Evaluation functions

			//- Return neighbour coupled given internal cell data
			tmp<Field<Type>> patchNeighbourField() const;

			//- Write
			virtual void write(Ostream&) const;

};

template<>
tmp<Field<scalar>> phaseLagCyclicFvPatchField<scalar>::phaseLaggedField() const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "phaseLagCyclicFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
