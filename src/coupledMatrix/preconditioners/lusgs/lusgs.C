/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa

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
