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

#include "phaseLagAMIFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makePatchFields(phaseLagAMI);

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::phaseLagAMIFvPatchField<scalar>::updateInterfaceMatrix
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
        this->phaseLagAMIPatch().phaseLagAMIPatch().neighbPatch().faceCells();

    solveScalarField pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    this->transformCoupleField(pnf, cmpt);

    const word& fieldName = this->internalField().name();

	if (phaseLagAMIPatch_.applyLowWeightCorrection())
	{
		solveScalarField pif(psiInternal, phaseLagAMIPatch_.faceCells());

	    if
	    (
	        reinterpret_cast<const void*>(&psiInternal)
	     == reinterpret_cast<const void*>(&this->primitiveField())
	    )
		{
			pnf = phaseLagAMIPatch_.interpolate(pnf, fieldName, pif);
		}
		else
		{
			pnf = phaseLagAMIPatch_.interpolate(pnf, pif);
		}
	}
	else
	{
	    if
	    (
	        reinterpret_cast<const void*>(&psiInternal)
	     == reinterpret_cast<const void*>(&this->primitiveField())
	    )
		{
			pnf = phaseLagAMIPatch_.interpolate(pnf, fieldName);
		}
		else
		{
			pnf = phaseLagAMIPatch_.interpolate(pnf);
		}
	}

	// Multiply the field by coefficients and add into the result
	this->addToInternalField(result, !add, coeffs, pnf);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
