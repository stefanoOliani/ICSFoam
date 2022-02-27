/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "jumpCyclicAMIFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makePatchFieldTypeNames(jumpCyclicAMI);

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::jumpCyclicAMIFvPatchField<scalar>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    const labelUList& nbrFaceCells =
        this->cyclicAMIPatch().cyclicAMIPatch().neighbPatch().faceCells();

    solveScalarField pnf(psiInternal, nbrFaceCells);

    pnf = this->cyclicAMIPatch().interpolate(pnf);

    // only apply jump to original field
    if
    (
        reinterpret_cast<const void*>(&psiInternal)
     == reinterpret_cast<const void*>(&this->primitiveField())
    )
    {
        Field<scalar> jf(this->jump());

        if (!this->cyclicAMIPatch().owner())
        {
            jf *= -1.0;
        }

        //pnf -= jf;
        forAll(pnf, i)
        {
            pnf[i] -= jf[i];
        }
    }

    // Transform according to the transformation tensors
    this->transformCoupleField(pnf, cmpt);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, coeffs, pnf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
