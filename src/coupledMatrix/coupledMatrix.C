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

#include "coupledMatrix.H"
#include "emptyFvPatch.H"
#include "transformFvPatchFields.H"
#include "blockFvOperators.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

	defineTypeNameAndDebug(coupledMatrix, 0);

	const label coupledMatrix::solver::defaultMaxIter_ = 1000;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * * * //

coupledMatrix::coupledMatrix
(
	const fvMesh& mesh,
	int nScalars,
	int nVectors,
	bool deltaForm
)
:
	nScal_(nScalars),
	nVect_(nVectors),
	deltaForm_(deltaForm),
    dSByS_(nScal_ * nScal_),
    dSByV_(nScal_ * nVect_),
    dVByS_(nVect_ * nScal_),
    dVByV_(nVect_ * nVect_),
    mesh_(mesh)
{}


// * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * * //

void coupledMatrix::matrixMul
(
    PtrList<volScalarField>& sVec, PtrList<volVectorField>& vVec,
    PtrList<scalarField>& sResult, PtrList<vectorField>& vResult
) const
{
    forN(nScal_,i) sResult[i] = Zero;
    forN(nVect_,i) vResult[i] = Zero;

    scalarField sTmp(mesh_.nCells());
    vectorField vTmp(mesh_.nCells());

    forN(nScal_, i)
    {
        forN(nScal_, j)
        {
            if (dSBySExists(i,j))
            {
                dSByS(i,j).Amul(sTmp, sVec[j], mesh_);
                sResult[i] += sTmp;
            }
        }
    }
    forN(nScal_, i)
    {
        forN(nVect_, j)
        {
            if (dSByVExists(i,j))
            {
                dSByV(i,j).Amul(sTmp, vVec[j], mesh_);
                sResult[i] += sTmp;
            }
        }
    }
    forN(nVect_, i)
    {
        forN(nScal_, j)
        {
            if (dVBySExists(i,j))
            {
                dVByS(i,j).Amul(vTmp, sVec[j], mesh_);
                vResult[i] += vTmp;
            }
        }
    }
    forN(nVect_, i)
    {
        forN(nVect_, j)
        {
            if (dVByVExists(i,j))
            {
                dVByV(i,j).Amul(vTmp, vVec[j], mesh_);
                vResult[i] += vTmp;
            }
        }
    }

}


void coupledMatrix::matrixMulNoDiag
(
    PtrList<volScalarField>& sVec, PtrList<volVectorField>& vVec,
    PtrList<scalarField>& sResult, PtrList<vectorField>& vResult
) const
{
    forN(nScal_,i) sResult[i] = Zero;
    forN(nVect_,i) vResult[i] = Zero;

    scalarField sTmp(mesh_.nCells(), 0.0);
    vectorField vTmp(mesh_.nCells(), Zero);

    forN(nScal_, i)
    {
        forN(nScal_, j)
        {
            if (dSBySExists(i,j))
            {
                dSByS(i,j).AmulNoDiag(sTmp, sVec[j], mesh_);
                sResult[i] += sTmp;
            }
        }
    }
    forN(nScal_, i)
    {
        forN(nVect_, j)
        {
            if (dSByVExists(i,j))
            {
                dSByV(i,j).AmulNoDiag(sTmp, vVec[j], mesh_);
                sResult[i] += sTmp;
            }
        }
    }
    forN(nVect_, i)
    {
        forN(nScal_, j)
        {
            if (dVBySExists(i,j))
            {
                dVByS(i,j).AmulNoDiag(vTmp, sVec[j], mesh_);
                vResult[i] += vTmp;
            }
        }
    }
    forN(nVect_, i)
    {
        forN(nVect_, j)
        {
            if (dVByVExists(i,j))
            {
                dVByV(i,j).AmulNoDiag(vTmp, vVec[j], mesh_);
                vResult[i] += vTmp;
            }
        }
    }
}


residualsIO coupledMatrix::solve
(
	const dictionary& solverControls,
	PtrList<volScalarField>& sW,
	PtrList<volVectorField>& vW
)
{
    PtrList<scalarField> scalarSources(nScal_);
    PtrList<vectorField> vectorSources(nVect_);

    forN(nScal_,i)
    {
    	scalarSources.set(i, new scalarField(sW[i].size(), 0.0));

    	forN(nScal_,j)
    	{
            scalarSources[i] += dSByS(i,j).source();
    	}

    	forN(nVect_,j)
    	{
    		scalarSources[i] += dSByV(i,j).source();
    	}
    }

    forN(nVect_,i)
    {
    	vectorSources.set(i, new vectorField(vW[i].size(), Zero));

    	forN(nScal_,j)
    	{
            vectorSources[i] += dVByS(i,j).source();
    	}

    	forN(nVect_,j)
    	{
    		vectorSources[i] += dVByV(i,j).source();
    	}
    }

    // Create and call solver
    autoPtr<coupledMatrix::solver> sol =
    	coupledMatrix::solver::New
        (
            solverControls,
            *this
        );

    Info<< sol().type()<< " : Solving for ( ";

	forAll(sW, j)
	{
	    Info << " ";
		Info << sW[j].name();
	}
	forAll(vW, j)
	{
	    Info << " ";
		Info << vW[j].name();
	}

    Info << " ) " << endl;

    residualsIO solverPerf(nScal_, nVect_);

    if (deltaForm_)
    {
    	solverPerf = sol->solveDelta(sW, vW, scalarSources, vectorSources);
    }
    else
    {
        solverPerf = sol->solve(sW, vW, scalarSources, vectorSources);
    }

    //Print performance
    solverPerf.print();

    forN(nScal_,i)
    {
        sW[i].correctBoundaryConditions();
    }
    forN(nVect_,i)
    {
        // Zero the variable in non-solved directions
        vector::labelType validComponents = this->mesh().solutionD(); //-1 for empty directions
        forAll(validComponents, cmpt)
        {
            if (validComponents[cmpt] == -1)
            {
                vW[i].replace(cmpt, dimensionedScalar("0", vW[i].dimensions(), 0.0));
            }
        }

        vW[i].correctBoundaryConditions();
    }

    return solverPerf;
}


residualsIO coupledMatrix::solve
(
	PtrList<volScalarField>& sW,
	PtrList<volVectorField>& vW
)
{
    const dictionary& dict = mesh().solutionDict().subDict("flowSolver");

    return solve(dict, sW, vW);

}

residualsIO coupledMatrix::solveForIncr
(
	PtrList<volScalarField>& sW,
	PtrList<volVectorField>& vW,
	PtrList<volScalarField>& dsW,
	PtrList<volVectorField>& dvW
)
{

	const dictionary& dict = mesh().solutionDict().subDict("flowSolver");

    PtrList<scalarField> scalarSources(nScal_);
    PtrList<vectorField> vectorSources(nVect_);

    forN(nScal_,i)
    {
    	scalarSources.set(i, new scalarField(sW[i].size(), 0.0));

    	forN(nScal_,j)
    	{
            scalarSources[i] += dSByS(i,j).source();
    	}

    	forN(nVect_,j)
    	{
    		scalarSources[i] += dSByV(i,j).source();
    	}
    }

    forN(nVect_,i)
    {
    	vectorSources.set(i, new vectorField(vW[i].size(), Zero));

    	forN(nScal_,j)
    	{
            vectorSources[i] += dVByS(i,j).source();
    	}

    	forN(nVect_,j)
    	{
    		vectorSources[i] += dVByV(i,j).source();
    	}
    }

    // Create and call solver
    autoPtr<coupledMatrix::solver> sol =
    	coupledMatrix::solver::New
        (
            dict,
            *this
        );

    Info<< sol().type()<< " : Solving for ( ";

	forAll(sW, j)
	{
	    Info << " ";
		Info << dsW[j].name();
	}
	forAll(vW, j)
	{
	    Info << " ";
		Info << dvW[j].name();
	}

    Info << " ) " << endl;

    residualsIO solverPerf(nScal_, nVect_);

    solverPerf = sol->solveDelta(sW, vW, scalarSources, vectorSources, dsW, dvW);

    //Print performance
    solverPerf.print();

    forN(nVect_,i)
    {
        // Zero the variable in non-solved directions
        vector::labelType validComponents = this->mesh().solutionD(); //-1 for empty directions
        forAll(validComponents, cmpt)
        {
            if (validComponents[cmpt] == -1)
            {
                dvW[i].replace(cmpt, dimensionedScalar("0", vW[i].dimensions(), 0.0));
            }
        }
    }

    return solverPerf;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
