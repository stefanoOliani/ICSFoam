/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

    Hrvoje Jasak, Wikki Ltd.  All rights reserved
    Fethi Tekin, All rights reserved.
    Oliver Borm, All rights reserved.

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

\*---------------------------------------------------------------------------*/

#include "overlapAMIFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makePatchFields(overlapAMI);

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::overlapAMIFvPatchField<scalar>::updateInterfaceMatrix
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
	const labelList& nbrFaceCells =
					lduAddr.patchAddr
					(
						overlapAMIPatch_.overlapAMIPatch().neighbPatchID()
					);

    solveScalarField pnf(psiInternal, nbrFaceCells);


	if (overlapAMIPatch_.applyLowWeightCorrection())
	{
		solveScalarField pif(psiInternal, overlapAMIPatch_.faceCells());

		pnf = overlapAMIPatch_.interpolate(pnf, pif);
	}
	else
	{
		pnf = overlapAMIPatch_.interpolate(pnf);
	}

	const labelUList& faceCells = lduAddr.patchAddr(patchId);

	// Multiply the field by coefficients and add into the result
	this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
