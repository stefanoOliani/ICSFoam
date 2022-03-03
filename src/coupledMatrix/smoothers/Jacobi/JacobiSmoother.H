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

#ifndef JacobiSmoother_H
#define JacobiSmoother_H

#include "coupledMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class JacobiSmoother Declaration
\*---------------------------------------------------------------------------*/

class JacobiSmoother
:
    public coupledMatrix::smoother
{
	//Store inverse of diagonal blocks
	Field<scalarSquareMatrix> InvDiag_;

	void calcInvDiag();


public:

    //- Runtime type information
    TypeName("Jacobi");


    // Constructors

        //- Construct from components
        JacobiSmoother
        (
            const coupledMatrix& matrix
        );


    // Member Functions

        //- Smooth for the given number of sweeps
        virtual void smooth
		(
			PtrList<volScalarField>& sW,
			PtrList<volVectorField>& vW,
			const PtrList<scalarField>& sSource,
			const PtrList<vectorField>& vSource,
			const label nSweeps
		) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
