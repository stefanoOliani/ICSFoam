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

#include "Jacobi.H"
#include "JacobiSmoother.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

	defineTypeNameAndDebug(Jacobi, 0);

    coupledMatrix::preconditioner::
        adddictionaryConstructorToTable<Jacobi>
        addJacobiDictionaryConstructorToTable_;

// * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * * * //

Jacobi::Jacobi
(
    const coupledMatrix::solver& sol,
    const dictionary&
)
:
coupledMatrix::preconditioner(sol)
{}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

void Jacobi::precondition
(
    PtrList<scalarField>& sVec,
    PtrList<vectorField>& vVec
) const
{
	const coupledMatrix& cMatrix = this->solver_.matrix();
	const int nScalar = cMatrix.nScal();
	const int nVector = cMatrix.nVect();

	// Allocate temp storage
    PtrList<volScalarField> dsW(nScalar);
    PtrList<volVectorField> dvW(nVector);
    forN(nScalar,i)
    {
        dsW.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "dScalar" + Foam::name(i),
					cMatrix.mesh().time().timeName(),
					cMatrix.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
				cMatrix.mesh(),
				dimensionedScalar
				(
					"dScalar",
					dimless,
					Zero
				)
            )
        );
    }
    forN(nVector,i)
    {
        dvW.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "dVector" + Foam::name(i),
					cMatrix.mesh().time().timeName(),
					cMatrix.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
				cMatrix.mesh(),
				dimensionedVector
				(
					"dVector",
					dimless,
					Zero
				)
            )
        );
    }

	JacobiSmoother JacSmoother(cMatrix);

	JacSmoother.smooth(dsW, dvW, sVec, vVec, 1);

    forN(nScalar,i)
	{
		sVec[i] = dsW[i].primitiveFieldRef();
	}
	forN(nVector,i)
	{
		vVec[i] = dvW[i].primitiveFieldRef();
	}

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
