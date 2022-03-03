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

#include "lusgs.H"
#include "emptyFvPatch.H"
#include "transformFvPatchFields.H"
#include "blockFvOperators.H"
#include "diagTensor.H"
#include "LUscalarMatrix.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

	defineTypeNameAndDebug(lusgs, 0);

    coupledMatrix::preconditioner::
        adddictionaryConstructorToTable<lusgs>
        addlusgsDictionaryConstructorToTable_;

// * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * * * //

lusgs::lusgs
(
    const coupledMatrix::solver& sol,
    const dictionary&
)
:
coupledMatrix::preconditioner(sol)
{
    // Generate a scalar diagonal coefficient based on the max of the diagonal
    // of the matrix

    rDiagCoeff_.set(new scalarField(sol.matrix().mesh().nCells(), GREAT));

    const coupledMatrix& cMatrix = sol.matrix();
    const int nScalar = cMatrix.nScal();
    const int nVector = cMatrix.nVect();

    forN(cMatrix.mesh().nCells(), celli)
    {
        forN(nScalar, i)
        {
            if (cMatrix.dSBySExists(i,i))
            {
                rDiagCoeff_()[celli] =
                    1.0/max
                    (
                        1.0/rDiagCoeff_()[celli],
                        mag(cMatrix.dSByS(i,i).diag()[celli])
                    );
            }
            else
            {
                FatalErrorInFunction
                    << "Diagonal S" << i << " of coupledMatrix not populated."
                    << exit(FatalError);
            }
        }
        forN(nVector, i)
        {
            if (cMatrix.dVByVExists(i,i))
            {
                const tensor& diag = cMatrix.dVByV(i,i).diag()[celli];

                rDiagCoeff_()[celli] =
                    1.0/max
                    (
                        1.0/rDiagCoeff_()[celli],
                        mag(diag.xx())
                    );
                rDiagCoeff_()[celli] =
                    1.0/max
                    (
                        1.0/rDiagCoeff_()[celli],
                        mag(diag.yy())
                    );
                rDiagCoeff_()[celli] =
                    1.0/max
                    (
                        1.0/rDiagCoeff_()[celli],
                        mag(diag.zz())
                    );
            }
            else
            {
                FatalErrorInFunction
                    << "Diagonal V" << i << " of coupledMatrix not populated."
                    << exit(FatalError);
            }
        }
        if (rDiagCoeff_()[celli] < VSMALL)
        {
            FatalErrorInFunction << "All diagonals of coupledMatrix are zero." << endl
                << exit(FatalError);
        }
    }
}


// * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * * * //

template <class Type1, class Type2>
void lusgs::forwardSweep
(
    const label celli,
    Field<Type1>& result,
	const blockFvMatrix<Type1, typename outerProduct<Type1, Type2>::type>& matrix,
    const Type2& delta
) const
{
    if (matrix.hasLower() || matrix.hasUpper())
    {
        const Field<typename outerProduct<Type1, Type2>::type>& lower =
            matrix.lower();
        const fvMesh& mesh = this->solver_.matrix().mesh();
        const labelUList& nei = mesh.neighbour();
        const cell& faces = mesh.cells()[celli];
        forAll(faces, facei)
        {
            const label faceI = faces[facei];
            if (faceI < mesh.nInternalFaces())
            {
                const label cellj = nei[faceI];
                if (cellj > celli)
                {
                    result[cellj] -= dot(lower[faceI], delta);
                }
            }
        }
    }
}


template <class Type1, class Type2>
void lusgs::reverseSweep
(
    const label celli,
    Field<Type1>& result,
	const blockFvMatrix<Type1, typename outerProduct<Type1, Type2>::type>& matrix,
    const Type2& delta
) const
{
    if (matrix.hasLower() || matrix.hasUpper())
    {
        const Field<typename outerProduct<Type1, Type2>::type>& upper =
            matrix.upper();
        const fvMesh& mesh = this->solver_.matrix().mesh();
        const labelUList& own = mesh.owner();
        const cell& faces = mesh.cells()[celli];
        forAll(faces, facei)
        {
            const label faceI = faces[facei];
            if (faceI < mesh.nInternalFaces())
            {
                const label cellj = own[faceI];
                if (cellj < celli)
                {
                    result[cellj] -= dot(upper[faceI], delta);
                }
            }
        }
    }
}


void lusgs::divideByDiagonal
(
    const label& celli,
    scalarList& dS,
    List<vector>& dV,
    const PtrList<scalarField>& sVec,
    const PtrList<vectorField>& vVec
) const
{
    const coupledMatrix& cMatrix = this->solver_.matrix();
    const int nScalar = cMatrix.nScal();
    const int nVector = cMatrix.nVect();

//    scalarSquareMatrix Jacobian_0(nScalar+3*nVector, 0.0);
//
//    const blockFvMatrix<scalar, scalar>& dContByRho = cMatrix.dSByS(0,0);
//    const blockFvMatrix<scalar, vector>& dContByRhoU = cMatrix.dSByV(0,0);
//    const blockFvMatrix<scalar, scalar>& dContByRhoE = cMatrix.dSByS(0,1);
//    const blockFvMatrix<vector, vector>& dMomByRho = cMatrix.dVByS(0,0);
//    const blockFvMatrix<vector, tensor>& dMomByRhoU = cMatrix.dVByV(0,0);
//    const blockFvMatrix<vector, vector>& dMomByRhoE = cMatrix.dVByS(0,1);
//    const blockFvMatrix<scalar, scalar>& dEnergyByRho = cMatrix.dSByS(1,0);
//    const blockFvMatrix<scalar, vector>& dEnergyByRhoU = cMatrix.dSByV(1,0);
//    const blockFvMatrix<scalar, scalar>& dEnergyByRhoE = cMatrix.dSByS(1,1);
//
//    scalar diagContByRho = dContByRho.diag()[celli];
//    vector diagContByRhoU = dContByRhoU.diag()[celli];
//    scalar diagContByRhoE = dContByRhoE.diag()[celli];
//
//    vector diagMomByRho = dMomByRho.diag()[celli];
//    tensor diagMomByRhoU = dMomByRhoU.diag()[celli];
//    vector diagMomByRhoE = dMomByRhoE.diag()[celli];
//
//    scalar diagEnergyByRho = dEnergyByRho.diag()[celli];
//    vector diagEnergyByRhoU = dEnergyByRhoU.diag()[celli];
//    scalar diagEnergyByRhoE = dEnergyByRhoE.diag()[celli];
//
//    Jacobian_0[0][0] = diagContByRho;
//    Jacobian_0[0][1] = diagContByRhoU.x();
//    Jacobian_0[0][2] = diagContByRhoU.y();
//    Jacobian_0[0][3] = diagContByRhoU.z();
//    Jacobian_0[0][4] = diagContByRhoE;
//
//    Jacobian_0[1][0] = diagMomByRho.x();
//    Jacobian_0[1][1] = diagMomByRhoU.xx();
//    Jacobian_0[1][2] = diagMomByRhoU.xy();
//    Jacobian_0[1][3] = diagMomByRhoU.xz();
//    Jacobian_0[1][4] = diagMomByRhoE.x();
//
//    Jacobian_0[2][0] = diagMomByRho.y();
//    Jacobian_0[2][1] = diagMomByRhoU.yx();
//    Jacobian_0[2][2] = diagMomByRhoU.yy();
//    Jacobian_0[2][3] = diagMomByRhoU.yz();
//    Jacobian_0[2][4] = diagMomByRhoE.y();
//
//    Jacobian_0[3][0] = diagMomByRho.z();
//    Jacobian_0[3][1] = diagMomByRhoU.zx();
//    Jacobian_0[3][2] = diagMomByRhoU.zy();
//    Jacobian_0[3][3] = diagMomByRhoU.zz();
//    Jacobian_0[3][4] = diagMomByRhoE.z();
//
//    Jacobian_0[4][0] = diagEnergyByRho;
//    Jacobian_0[4][1] = diagEnergyByRhoU.x();
//    Jacobian_0[4][2] = diagEnergyByRhoU.y();
//    Jacobian_0[4][3] = diagEnergyByRhoU.z();
//    Jacobian_0[4][4] = diagEnergyByRhoE;
//
//    LUscalarMatrix Jacobian(Jacobian_0);
//
//    scalarSquareMatrix invJac(nScalar+3*nVector,0.0);
//    Jacobian.inv(invJac);
//
//    scalarList variables(nScalar+3*nVector, 0.0);
//    scalarList result(nScalar+3*nVector, 0.0);
//
//    variables[0] = sVec[0][celli];
//    variables[1] = vVec[0][celli].x();
//    variables[2] = vVec[0][celli].y();
//    variables[3] = vVec[0][celli].z();
//    variables[4] = sVec[1][celli];
//
//    for (int i=0; i<5; i++)
//    {
//    	for (int j=0; j<5; j++)
//    	{
//    		result[i] += invJac[i][j]*variables[j];
//    	}
//    }
//
//    dS[0] = result[0];
//    dV[0].x() = result[1];
//    dV[0].y() = result[2];
//    dV[0].z() = result[3];
//    dS[1] = result[4];


//    scalarSquareMatrix Jacobian_0(nScalar+3*nVector, 0.0);
//
//    const blockFvMatrix<scalar, scalar>& dContByP = cMatrix.dSByS(0,0);
//    const blockFvMatrix<scalar, vector>& dContByU = cMatrix.dSByV(0,0);
//    const blockFvMatrix<vector, vector>& dMomByP = cMatrix.dVByS(0,0);
//    const blockFvMatrix<vector, tensor>& dMomByU = cMatrix.dVByV(0,0);
//
//    scalar diagContByP = dContByP.diag()[celli];
//    vector diagContByU = dContByU.diag()[celli];
//
//    vector diagMomByP = dMomByP.diag()[celli];
//    tensor diagMomByU = dMomByU.diag()[celli];
//
//    Jacobian_0[0][0] = diagContByP;
//    Jacobian_0[0][1] = diagContByU.x();
//    Jacobian_0[0][2] = diagContByU.y();
//    Jacobian_0[0][3] = diagContByU.z();
//
//    Jacobian_0[1][0] = diagMomByP.x();
//    Jacobian_0[1][1] = diagMomByU.xx();
//    Jacobian_0[1][2] = diagMomByU.xy();
//    Jacobian_0[1][3] = diagMomByU.xz();
//
//    Jacobian_0[2][0] = diagMomByP.y();
//    Jacobian_0[2][1] = diagMomByU.yx();
//    Jacobian_0[2][2] = diagMomByU.yy();
//    Jacobian_0[2][3] = diagMomByU.yz();
//
//    Jacobian_0[3][0] = diagMomByP.z();
//    Jacobian_0[3][1] = diagMomByU.zx();
//    Jacobian_0[3][2] = diagMomByU.zy();
//    Jacobian_0[3][3] = diagMomByU.zz();
//
//    LUscalarMatrix Jacobian(Jacobian_0);
//
//    scalarSquareMatrix invJac(nScalar+3*nVector,0.0);
//    Jacobian.inv(invJac);
//
//    scalarList variables(nScalar+3*nVector, 0.0);
//    scalarList result(nScalar+3*nVector, 0.0);
//
//    variables[0] = sVec[0][celli];
//    variables[1] = vVec[0][celli].x();
//    variables[2] = vVec[0][celli].y();
//    variables[3] = vVec[0][celli].z();
//
//    for (int i=0; i<4; i++)
//    {
//    	for (int j=0; j<4; j++)
//    	{
//    		result[i] += invJac[i][j]*variables[j];
//    	}
//    }
//
//    dS[0] = result[0];
//    dV[0].x() = result[1];
//    dV[0].y() = result[2];
//    dV[0].z() = result[3];


    forN(nScalar, i)
    {
        dS[i] = rDiagCoeff_()[celli] * sVec[i][celli];
    }
    forN(nVector, i)
    {
        dV[i] = rDiagCoeff_()[celli] * vVec[i][celli];
    }

}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

void lusgs::precondition
(
    PtrList<scalarField>& sVec,
    PtrList<vectorField>& vVec
) const
{
    const coupledMatrix& cMatrix = this->solver_.matrix();
    const int nScalar = cMatrix.nScal();
    const int nVector = cMatrix.nVect();

    const fvMesh& mesh = this->solver_.matrix().mesh();

//    // Residual is still in strong form
//    forN(nScalar, i) sVec[i] *= mesh.V();
//    forN(nVector, i) vVec[i] *= mesh.V();

    // Lower sweep: D \Delta W* = ( R - L \Delta W* )
    forAll(mesh.cells(), celli)
    {
        scalarList dSStar(nScalar);
        List<vector> dVStar(nVector);
        divideByDiagonal(celli, dSStar, dVStar, sVec, vVec);

        // Distribute to future cells
        forN(nScalar, i)
        {
            forN(nScalar, j)
            {
                if (cMatrix.dSBySExists(i,j))
                {
                    forwardSweep
                    (
                        celli,
                        sVec[i],
						cMatrix.dSByS(i,j),
                        dSStar[j]
                    );
                }
            }
        }
        forN(nScalar, i)
        {
            forN(nVector, j)
            {
                if (cMatrix.dSByVExists(i,j))
                {
                    forwardSweep
                    (
                        celli,
                        sVec[i],
						cMatrix.dSByV(i,j),
                        dVStar[j]
                    );
                }
            }
        }
        forN(nVector, i)
        {
            forN(nScalar, j)
            {
                if (cMatrix.dVBySExists(i,j))
                {
                    forwardSweep
                    (
                        celli,
                        vVec[i],
						cMatrix.dVByS(i,j),
                        dSStar[j]
                    );
                }
            }
        }
        forN(nVector, i)
        {
            forN(nVector, j)
            {
                if (cMatrix.dVByVExists(i,j))
                {
                    forwardSweep
                    (
                        celli,
                        vVec[i],
						cMatrix.dVByV(i,j),
                        dVStar[j]
                    );
                }
            }
        }
    }

    // Upper sweep: \Delta W = rD( D \Delta W* - U \Delta W )
    forAllReverse(mesh.cells(), celli)
    {
        scalarList dS(nScalar);
        List<vector> dV(nVector);
        divideByDiagonal(celli, dS, dV, sVec, vVec);
        forAll(sVec,i) sVec[i][celli] = dS[i];
        forAll(vVec,i) vVec[i][celli] = dV[i];

        // Distribute to future cells
        // For symmetric matrices, upper() returns lower
        forN(nScalar, i)
        {
            forN(nScalar, j)
            {
                if (cMatrix.dSBySExists(i,j))
                {
                    reverseSweep
                    (
                        celli,
                        sVec[i],
						cMatrix.dSByS(i,j),
                        dS[j]
                    );
                }
            }
        }
        forN(nScalar, i)
        {
            forN(nVector, j)
            {
                if (cMatrix.dSByVExists(i,j))
                {
                    reverseSweep
                    (
                        celli,
                        sVec[i],
						cMatrix.dSByV(i,j),
                        dV[j]
                    );
                }
            }
        }
        forN(nVector, i)
        {
            forN(nScalar, j)
            {
                if (cMatrix.dVBySExists(i,j))
                {
                    reverseSweep
                    (
                        celli,
                        vVec[i],
						cMatrix.dVByS(i,j),
                        dS[j]
                    );
                }
            }
        }
        forN(nVector, i)
        {
            forN(nVector, j)
            {
                if (cMatrix.dVByVExists(i,j))
                {
                    reverseSweep
                    (
                        celli,
                        vVec[i],
						cMatrix.dVByV(i,j),
                        dV[j]
                    );
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
