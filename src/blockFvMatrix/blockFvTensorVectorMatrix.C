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

#include "blockFvTensorVectorMatrix.H"
#include "fvScalarMatrix.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::blockFvMatrix<Foam::vector,Foam::tensor>::insertEquation
(
	fvVectorMatrix& matrix
)
{
// Save a copy for different components
    scalarField& diag = matrix.diag();
    scalarField saveDiag(diag);

    // Add source boundary contribution
    vectorField& source = matrix.source();

    addBoundarySource(matrix, source, false);

    // Get reference to this source field of block system
    vectorField& b = this->source();

    tensorField& blockDiag = this->diag();

    for (direction cmptI = 0; cmptI < 3; cmptI++)
    {
        addBoundaryDiag(matrix, diag, cmptI);
        scalarField sourceCmpt(source.component(cmptI));

        forAll (diag, cellI)
        {
		blockDiag[cellI](cmptI, cmptI) = diag[cellI];
		b[cellI].component(cmptI) += sourceCmpt[cellI];
        }

        // Reset diagonal
         diag = saveDiag;
    }

    if (matrix.hasUpper())
    {
        const scalarField& upper = matrix.upper();

	tensorField& blockUpper = this->upper();

        for (direction cmptI = 0; cmptI < 3; cmptI++)
        {
            forAll (upper, faceI)
            {
                blockUpper[faceI](cmptI, cmptI) = upper[faceI];
            }
        }
    }

    if (matrix.hasLower())
    {
        const scalarField& lower = matrix.lower();

	tensorField& blockLower = this->lower();

        for (direction cmptI = 0; cmptI < 3; cmptI++)
        {
            forAll (lower, faceI)
            {
		blockLower[faceI](cmptI, cmptI) = lower[faceI];
            }
        }
    }

    const volVectorField& psi = matrix.psi();

    const fvMesh& mesh = psi.mesh();

    this->interfacesUpper().resize(mesh.boundary().size());
    this->interfacesLower().resize(mesh.boundary().size());

    forAll (psi.boundaryField(), patchI)
    {
	const fvPatchField<vector>& pf = psi.boundaryField()[patchI];

	const vectorField& icp = matrix.internalCoeffs()[patchI];
	const vectorField& bcp = matrix.boundaryCoeffs()[patchI];

        this->interfacesUpper().set(patchI, new tensorField(icp.size(), tensor(Zero)));
        this->interfacesLower().set(patchI, new tensorField(icp.size(), tensor(Zero)));

	tensorField& pcoupleUpper = this->interfacesUpper()[patchI];
	tensorField& pcoupleLower = this->interfacesLower()[patchI];

	for (direction cmptI = 0; cmptI < 3; cmptI++)
	{
		scalarField icpCmpt = icp.component(cmptI);
		scalarField bcpCmpt = bcp.component(cmptI);

		forAll (pf, faceI)
		{
			pcoupleUpper[faceI](cmptI, cmptI) = -bcpCmpt[faceI];
			pcoupleLower[faceI](cmptI, cmptI) = icpCmpt[faceI];
		}
	}
    }
}

// ************************************************************************* //
