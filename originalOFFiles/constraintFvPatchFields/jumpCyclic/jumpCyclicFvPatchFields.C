/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include "jumpCyclicFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makePatchFieldTypeNames(jumpCyclic);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::jumpCyclicFvPatchField<Foam::scalar>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    solveScalarField pnf(this->size());

    const labelUList& nbrFaceCells =
        lduAddr.patchAddr
        (
            this->cyclicPatch().neighbPatchID()
        );

    // only apply jump to original field
    if
    (
        reinterpret_cast<const void*>(&psiInternal)
     == reinterpret_cast<const void*>(&this->primitiveField())
    )
    {
        Field<scalar> jf(this->jump());

        if (!this->cyclicPatch().owner())
        {
            jf *= -1.0;
        }

        forAll(*this, facei)
        {
            pnf[facei] = psiInternal[nbrFaceCells[facei]] - jf[facei];
        }
    }
    else
    {
        forAll(*this, facei)
        {
            pnf[facei] = psiInternal[nbrFaceCells[facei]];
        }
    }

    // Transform according to the transformation tensors
    this->transformCoupleField(pnf, cmpt);

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}


template<>
void Foam::jumpCyclicFvPatchField<Foam::vector>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    solveScalarField pnf(this->size());

    const labelUList& nbrFaceCells =
        lduAddr.patchAddr
        (
            this->cyclicPatch().neighbPatchID()
        );

    const Field<vector>& iField = this->primitiveField();

    // only apply jump to original field
    if
    (
        reinterpret_cast<const void*>(&psiInternal)
     == reinterpret_cast<const void*>(&(iField.component(cmpt).ref()))
    )
    {
        Field<vector> jf(this->jump());

        if (!this->cyclicPatch().owner())
        {
            jf *= -1.0;
        }

        forAll(*this, facei)
        {
            pnf[facei] =
                psiInternal[nbrFaceCells[facei]]
              - jf[facei].component(cmpt);
        }
    }
    else
    {
        forAll(*this, facei)
        {
            pnf[facei] = psiInternal[nbrFaceCells[facei]];
        }
    }

    // Transform according to the transformation tensors
    this->transformCoupleField(pnf, cmpt);

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}

// ************************************************************************* //
