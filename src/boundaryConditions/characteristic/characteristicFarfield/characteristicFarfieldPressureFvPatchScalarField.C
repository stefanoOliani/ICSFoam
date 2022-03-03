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

#include "characteristicFarfieldPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::characteristicFarfieldPressureFvPatchScalarField::
characteristicFarfieldPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    characteristicBase(p)
{
    refValue() = patchInternalField();
    refGrad() = 0;
    valueFraction() = 1;
}


Foam::characteristicFarfieldPressureFvPatchScalarField::
characteristicFarfieldPressureFvPatchScalarField
(
    const characteristicFarfieldPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    characteristicBase(ptf, p, mapper)
{}


Foam::characteristicFarfieldPressureFvPatchScalarField::
characteristicFarfieldPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    characteristicBase(p, dict)
{

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
    }

    refValue() = *this;
    refGrad() = 0;
    valueFraction() = 1;
}


Foam::characteristicFarfieldPressureFvPatchScalarField::
characteristicFarfieldPressureFvPatchScalarField
(
    const characteristicFarfieldPressureFvPatchScalarField& sfspvf
)
:
    mixedFvPatchScalarField(sfspvf),
    characteristicBase(sfspvf)
{}


Foam::characteristicFarfieldPressureFvPatchScalarField::
characteristicFarfieldPressureFvPatchScalarField
(
    const characteristicFarfieldPressureFvPatchScalarField& sfspvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(sfspvf, iF),
    characteristicBase(sfspvf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::characteristicFarfieldPressureFvPatchScalarField::updateCoeffs()
{
    if (!size() || updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const fluidThermo& thermo =
        mesh.lookupObject<fluidThermo>("thermophysicalProperties");

    tmp< volScalarField > gamma = thermo.gamma();
    const fvPatchField<scalar>& pgamma =
        gamma->boundaryField()[patch().index()];

    const fvPatchField<scalar>& ppsi =
        thermo.psi().boundaryField()[patch().index()];

    const fvPatchField<vector>& pU =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const fvsPatchField<scalar>& pphi =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);
    const fvPatchField<scalar>& prho =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    scalarField& pp = refValue();
    valueFraction() = 1;
    refGrad() = 0;

    // get the near patch internal cell values
    const scalarField p(patchInternalField());
    const vectorField U(pU.patchInternalField());

    // Patch outward pointing face unit vector (Same convention as Blazek)
    const vectorField pn(patch().nf());

    // Patch normal Mach number
    const scalarField pc(sqrt(pgamma/ppsi));
    const scalarField pM(pphi/(prho*patch().magSf()*pc));

    // Reference values (Blazek suggests using internal values at cell centres)
    const scalarField cO
    (
        sqrt(pgamma.patchInternalField()/ppsi.patchInternalField())
    );
    const scalarField rhoO(prho.patchInternalField());

    // Set the patch boundary condition based on the Mach number and direction
    // of the flow dictated by the boundary/free-stream pressure difference

    forAll(pp, facei)
    {
        if (pM[facei] <= -1.0)                       // Supersonic inflow
        {
            pp[facei] = pRef_;
        }
        else if (pM[facei] >= 1.0)                   // Supersonic outflow
        {
            valueFraction()[facei] = 0;
        }
        else if (pM[facei] <= 0.0)                   // Subsonic inflow
        {
            valueFraction()[facei] = 0.5;
            pp[facei] =
                pRef_
              - (rhoO[facei]*cO[facei]) * ((URef_ - U[facei]) & pn[facei]);
        }
        else                                         // Subsonic outflow
        {
            valueFraction()[facei] = 0.5;
            pp[facei] =
                pRef_
              - (rhoO[facei]*cO[facei]) * ((URef_ - U[facei]) & pn[facei]);
        }

    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::characteristicFarfieldPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    characteristicBase::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        characteristicFarfieldPressureFvPatchScalarField
    );
}

// ************************************************************************* //
