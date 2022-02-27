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

#ifndef characteristicPressureInletOutletPressureFvPatchScalarFields_H
#define characteristicPressureInletOutletPressureFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"
#include "characteristicBase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class characteristicPressureInletOutletPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class characteristicPressureInletOutletPressureFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public characteristicBase
{
public:

    //- Runtime type information
    TypeName("characteristicPressureInletOutletPressure");


    // Constructors

        //- Construct from patch and internal field
        characteristicPressureInletOutletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        characteristicPressureInletOutletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given characteristicPressureInletOutletPressureFvPatchScalarField
        //  onto a new patch
        characteristicPressureInletOutletPressureFvPatchScalarField
        (
            const characteristicPressureInletOutletPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        characteristicPressureInletOutletPressureFvPatchScalarField
        (
            const characteristicPressureInletOutletPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new characteristicPressureInletOutletPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        characteristicPressureInletOutletPressureFvPatchScalarField
        (
            const characteristicPressureInletOutletPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new characteristicPressureInletOutletPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
