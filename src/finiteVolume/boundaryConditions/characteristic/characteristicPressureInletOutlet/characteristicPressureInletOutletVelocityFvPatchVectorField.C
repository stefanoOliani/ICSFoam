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

#include "characteristicPressureInletOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::characteristicPressureInletOutletVelocityFvPatchVectorField::
characteristicPressureInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    characteristicBase(p)
{
    refValue() = patchInternalField();
    refGrad() = vector::zero;
    valueFraction() = 1;
}


Foam::characteristicPressureInletOutletVelocityFvPatchVectorField::
characteristicPressureInletOutletVelocityFvPatchVectorField
(
    const characteristicPressureInletOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    characteristicBase(ptf, p, mapper)
{}


Foam::characteristicPressureInletOutletVelocityFvPatchVectorField::
characteristicPressureInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    characteristicBase(p, dict)
{
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<vector>::operator=(patchInternalField());
    }

    refValue() = *this;
    refGrad() = vector::zero;
    valueFraction() = 1;
}


Foam::characteristicPressureInletOutletVelocityFvPatchVectorField::
characteristicPressureInletOutletVelocityFvPatchVectorField
(
    const characteristicPressureInletOutletVelocityFvPatchVectorField& sfspvf
)
:
    mixedFvPatchVectorField(sfspvf),
    characteristicBase(sfspvf)
{}


Foam::characteristicPressureInletOutletVelocityFvPatchVectorField::
characteristicPressureInletOutletVelocityFvPatchVectorField
(
    const characteristicPressureInletOutletVelocityFvPatchVectorField& sfspvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(sfspvf, iF),
    characteristicBase(sfspvf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::characteristicPressureInletOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (!size() || updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const fluidThermo& thermo =
        mesh.lookupObject<fluidThermo>("thermophysicalProperties");

    tmp< volScalarField > gamma = thermo.gamma();
    const fvPatchField<scalar>&  pgamma =
        gamma->boundaryField()[patch().index()];

    const fvPatchField<scalar>& ppsi =
        thermo.psi().boundaryField()[patch().index()];

    const fvPatchField<scalar>& pp =
        patch().lookupPatchField<volScalarField, scalar>(pName_);

    const fvsPatchField<scalar>& pphi =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);
    const fvPatchField<scalar>& prho =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    vectorField& Up = refValue();
    valueFraction() = 1;
    refGrad() = Zero;

    // get the near patch internal cell values
    const vectorField U(patchInternalField());
    const scalarField p(pp.patchInternalField());

    // Patch outward pointing unit vector (Same convention as Blazek)
    const vectorField pn(patch().nf());

    // Patch normal Mach number
    const scalarField pc(sqrt(pgamma/ppsi));
    const scalarField pM(pphi/(prho*patch().magSf()*pc));

    // Reference values (Blazek suggests using internal values at cell centres)
    scalarField cO
    (
        sqrt(pgamma.patchInternalField()/ppsi.patchInternalField())
    );
    scalarField rhoO(prho.patchInternalField());

    const scalarField& deltaCoeffs(patch().deltaCoeffs());

    // Set the patch boundary condition based on the Mach number and direction
    // of the flow dictated by the boundary/free-stream pressure difference

    forAll(Up, facei)
    {
        if (pM[facei] <= -1.0)                       // Supersonic inflow
        {
            Up[facei] = URef_;
        }
        else if (pM[facei] >= 1.0)                  // Supersonic outflow
        {
            valueFraction()[facei] = 0;
        }
        else if (pM[facei] <= 0.0)                   // Subsonic inflow
        {
            // Prescribed p and u_t
            // Extrapolated u_n
            valueFraction()[facei] = 0;
            vector Utang = U[facei]-(U[facei]&pn[facei])*pn[facei];
            vector Ureftang = URef_-(URef_&pn[facei])*pn[facei];
            refGrad()[facei] =
                (
                    -Utang + Ureftang
                )*deltaCoeffs[facei];
        }
        else                                         // Subsonic outflow
        {
            // Extrapolated U
            valueFraction()[facei] = 0;
            refGrad()[facei] = vector::zero;
        }

    }

    mixedFvPatchVectorField::updateCoeffs();
}


void Foam::characteristicPressureInletOutletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    characteristicBase::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        characteristicPressureInletOutletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
