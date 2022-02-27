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

#include "characteristicPressureInletOutletTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::characteristicPressureInletOutletTemperatureFvPatchScalarField::
characteristicPressureInletOutletTemperatureFvPatchScalarField
(
    const fvPatch& t,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(t, iF),
    characteristicBase(t)
{
    refValue() = patchInternalField();
    refGrad() = 0;
    valueFraction() = 1;
}


Foam::characteristicPressureInletOutletTemperatureFvPatchScalarField::
characteristicPressureInletOutletTemperatureFvPatchScalarField
(
    const characteristicPressureInletOutletTemperatureFvPatchScalarField& ptf,
    const fvPatch& t,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, t, iF, mapper),
    characteristicBase(ptf, t, mapper)
{}


Foam::characteristicPressureInletOutletTemperatureFvPatchScalarField::
characteristicPressureInletOutletTemperatureFvPatchScalarField
(
    const fvPatch& t,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(t, iF),
    characteristicBase(t, dict)
{

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, t.size())
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


Foam::characteristicPressureInletOutletTemperatureFvPatchScalarField::
characteristicPressureInletOutletTemperatureFvPatchScalarField
(
    const characteristicPressureInletOutletTemperatureFvPatchScalarField& sfspvf
)
:
    mixedFvPatchScalarField(sfspvf),
    characteristicBase(sfspvf)
{}


Foam::characteristicPressureInletOutletTemperatureFvPatchScalarField::
characteristicPressureInletOutletTemperatureFvPatchScalarField
(
    const characteristicPressureInletOutletTemperatureFvPatchScalarField& sfspvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(sfspvf, iF),
    characteristicBase(sfspvf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::characteristicPressureInletOutletTemperatureFvPatchScalarField::updateCoeffs()
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

    const fvPatchField<scalar>& pp =
        patch().lookupPatchField<volScalarField, scalar>(pName_);

    const fvPatchField<vector>& pU =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const fvsPatchField<scalar>& pphi =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);
    const fvPatchField<scalar>& prho =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    scalarField& pT = refValue();
    valueFraction() = 1;
    refGrad() = 0;

    // get the near patch internal cell values
    const scalarField T(patchInternalField());
    const scalarField p(pp.patchInternalField());
    const vectorField U(pU.patchInternalField());

    // Patch outward pointing unit vector (Same convention as Blazek)
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

    // Need effective R of the free-stream flow
    const scalarField Reff(p/(rhoO*patchInternalField()));
    const scalarField rhoInf(pRef_/(Reff*TRef_));

    // Set the patch boundary condition based on the Mach number and direction
    // of the flow dictated by the boundary/free-stream pressure difference

    forAll(pT, facei)
    {
        if (pM[facei] <= -1.0)                       // Supersonic inflow
        {
            pT[facei] = TRef_;
        }
        else if (pM[facei] >= 1.0)                   // Supersonic outflow
        {
            valueFraction()[facei] = 0;
        }
        else if (pM[facei] <= 0.0)                   // Subsonic inflow
        {
            pT[facei] = TRef_;
        }
        else                                         // Subsonic outflow
        {
            // Extrapolated T
            valueFraction()[facei] = 0;
        }
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::characteristicPressureInletOutletTemperatureFvPatchScalarField::write(Ostream& os) const
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
        characteristicPressureInletOutletTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
