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

#include "blockFvScalarScalarMatrix.H"
#include "fvScalarMatrix.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::blockFvMatrix<Foam::scalar,Foam::scalar>::insertEquation
(
	fvScalarMatrix& matrix
)
{
// Save a copy for different components
    scalarField& diag = matrix.diag();
    scalarField saveDiag(diag);

    // Add source boundary contribution
    scalarField& source = matrix.source();
    addBoundarySource(matrix,source, false);

    // Get reference to this source field of block system
    scalarField& b = this->source();

    scalarField& blockDiag = this->diag();

    addBoundaryDiag(matrix,diag,0);

    forAll (diag, cellI)
    {
    	blockDiag[cellI] = diag[cellI];
    	b[cellI] += source[cellI];
    }
    
    // Reset diagonal
    diag = saveDiag;

    if (matrix.hasUpper())
    {
        const scalarField& upper = matrix.upper();
	scalarField& blockUpper = this->upper();

	forAll (upper, faceI)
	{
		blockUpper[faceI] = upper[faceI];
	}
    }

    if (matrix.hasLower())
    {
        const scalarField& lower = matrix.lower();

	scalarField& blockLower = this->lower();

	forAll (lower, faceI)
	{
		blockLower[faceI] = lower[faceI];
	}
    }

    const volScalarField& psi = matrix.psi();

    const fvMesh& mesh = psi.mesh();

    this->interfacesUpper().resize(mesh.boundary().size());
    this->interfacesLower().resize(mesh.boundary().size());

    forAll (psi.boundaryField(), patchI)
    {
	const scalarField& icp = matrix.internalCoeffs()[patchI];
	const scalarField& bcp = matrix.boundaryCoeffs()[patchI];

        this->interfacesUpper().set(patchI, new scalarField(-bcp));
        this->interfacesLower().set(patchI, new scalarField(icp));
    }
}

// ************************************************************************* //
