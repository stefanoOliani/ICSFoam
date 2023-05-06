/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

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

#include "JacobiSmoother.H"
#include "diagTensor.H"
#include "LUscalarMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(JacobiSmoother, 0);

    coupledMatrix::smoother::adddictionaryConstructorToTable<JacobiSmoother>
        addJacobiSmootheDictionaryConstructorToTable_;
}


// * * * * * * * * *  *  * * Private member functions  * * * *  * * * * * * * //

void Foam::JacobiSmoother::calcInvDiag()
{
	const coupledMatrix& cMatrix = this->matrix();
	const int nScalar = this->matrix().nScal();
	const int nVector = this->matrix().nVect();
	const int nVariables = nScalar + 3*nVector;
	label nCells = this->matrix().mesh().nCells();

	for (label celli=0; celli<nCells; celli++)
	{
		scalarSquareMatrix Jacobian_0(nVariables, 0.0);

		forN(nScalar, s)
		{
			forN(nScalar, ns)
			{
				Jacobian_0[s][ns] = cMatrix.dSByS(s,ns).diag()[celli];
			}

			forN(nVector, nv)
			{
				Jacobian_0[s][nScalar + 3*nv] = cMatrix.dSByV(s,nv).diag()[celli].x();
				Jacobian_0[s][nScalar + 3*nv +1] = cMatrix.dSByV(s,nv).diag()[celli].y();
				Jacobian_0[s][nScalar + 3*nv +2] = cMatrix.dSByV(s,nv).diag()[celli].z();
			}
		}

		forN(nVector, v)
		{
			forN(nScalar, ns)
			{
				Jacobian_0[nScalar+3*v][ns] = cMatrix.dVByS(v,ns).diag()[celli].x();
				Jacobian_0[nScalar+3*v+1][ns] = cMatrix.dVByS(v,ns).diag()[celli].y();
				Jacobian_0[nScalar+3*v+2][ns] = cMatrix.dVByS(v,ns).diag()[celli].z();
			}

			forN(nVector, nv)
			{
				Jacobian_0[nScalar+3*v][nScalar + nv] = cMatrix.dVByV(v,nv).diag()[celli].xx();
				Jacobian_0[nScalar+3*v+1][nScalar + nv] = cMatrix.dVByV(v,nv).diag()[celli].yx();
				Jacobian_0[nScalar+3*v+2][nScalar + nv] = cMatrix.dVByV(v,nv).diag()[celli].zx();

				Jacobian_0[nScalar+3*v][nScalar + nv + 1] = cMatrix.dVByV(v,nv).diag()[celli].xy();
				Jacobian_0[nScalar+3*v+1][nScalar + nv + 1] = cMatrix.dVByV(v,nv).diag()[celli].yy();
				Jacobian_0[nScalar+3*v+2][nScalar + nv + 1] = cMatrix.dVByV(v,nv).diag()[celli].zy();

				Jacobian_0[nScalar+3*v][nScalar + nv +2] = cMatrix.dVByV(v,nv).diag()[celli].xz();
				Jacobian_0[nScalar+3*v+1][nScalar + nv + 2] = cMatrix.dVByV(v,nv).diag()[celli].yz();
				Jacobian_0[nScalar+3*v+2][nScalar + nv + 2] = cMatrix.dVByV(v,nv).diag()[celli].zz();
			}
		}

		LUscalarMatrix Jacobian(Jacobian_0);

		scalarSquareMatrix invJac(nVariables,0.0);
		Jacobian.inv(invJac);

		InvDiag_[celli] = invJac;
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::JacobiSmoother::JacobiSmoother
(
    const coupledMatrix& matrix
)
:
    coupledMatrix::smoother
    (
        matrix
    ),
	InvDiag_(matrix.mesh().cells().size())
{
	calcInvDiag();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::JacobiSmoother::smooth
(
	PtrList<volScalarField>& sW,
	PtrList<volVectorField>& vW,
	const PtrList<scalarField>& sSource,
	const PtrList<vectorField>& vSource,
	const label nSweeps
) const
{
	const int nScalar = this->matrix().nScal();
	const int nVector = this->matrix().nVect();
	const int nVariables = nScalar + 3*nVector;

    // Allocate temp storage
    PtrList<scalarField> sTmp(nScalar);
    PtrList<vectorField> vTmp(nVector);

    forN(nScalar,j)
    {
        sTmp.set(j, new scalarField(sW[j].size(), 0.0));
    }

    forN(nVector,j)
    {
    	vTmp.set(j, new vectorField(vW[j].size(), Zero));
    }

	label nCells = this->matrix().mesh().nCells();

    for (label sweep=0; sweep<nSweeps; sweep++)
    {
        // Calculate initial vector
        matrix_.matrixMulNoDiag(sW, vW, sTmp, vTmp);

        forN(nScalar,j)
        {
           sTmp[j] = -(sTmp[j] - sSource[j]);
        }

        forN(nVector,j)
        {
        	vTmp[j] = -(vTmp[j] - vSource[j]);
        }

        for (label celli=0; celli<nCells; celli++)
        {
        	scalarList variables(nVariables, 0.0);
        	scalarList result(nVariables, 0.0);

        	forN(nScalar, ns)
        	{
        		variables[ns] = sTmp[ns][celli];
        	}
        	forN(nVector, nv)
        	{
        		variables[nScalar + 3*nv] = vTmp[nv][celli].x();
        		variables[nScalar + 3*nv + 1] = vTmp[nv][celli].y();
        		variables[nScalar + 3*nv + 2] = vTmp[nv][celli].z();
        	}

        	for (int i=0; i<nVariables; i++)
        	{
        	  	for (int j=0; j<nVariables; j++)
        	 	{
        	   		result[i] += InvDiag_[celli][i][j]*variables[j];
        	   	}
        	}

        	forN(nScalar, ns)
        	{
        		sW[ns].primitiveFieldRef()[celli] = result[ns];
        	}
        	forN(nVector, nv)
        	{
        		vW[nv].primitiveFieldRef()[celli].x() = result[nScalar + 3*nv];
        		vW[nv].primitiveFieldRef()[celli].y() = result[nScalar + 3*nv + 1];
        		vW[nv].primitiveFieldRef()[celli].z() = result[nScalar + 3*nv + 2];
        	}
        }
    }
}

// ************************************************************************* //
