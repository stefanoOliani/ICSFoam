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


#include "smoothSolverCoupled.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

	defineTypeNameAndDebug(smoothSolverCoupled, 0);

	coupledMatrix::solver::adddictionaryConstructorToTable<smoothSolverCoupled>
    	addsmoothSolverCoupledDictionaryConstructorToTable_;

// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

smoothSolverCoupled::smoothSolverCoupled
(
    const dictionary& dict,
    const coupledMatrix& matrix
)
:
	coupledMatrix::solver(typeName, dict, matrix)
{
	readControls();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::smoothSolverCoupled::readControls()
{
    coupledMatrix::solver::readControls();
    nSweeps_ = controlDict_.getOrDefault<label>("nSweeps", 1);
}

residualsIO smoothSolverCoupled::solve
(
    PtrList<volScalarField>& sW, PtrList<volVectorField>& vW,
    const PtrList<scalarField>& sSource, const PtrList<vectorField>& vSource
) const
{
    const int nScalar = this->matrix().nScal();
    const int nVector = this->matrix().nVect();

	residualsIO solverPerf(nScalar,nVector);

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

    scalarList sNormFactor(nScalar);
    List<vector> vNormFactor(nVector);

    this->normFactors(sW,vW,sSource,vSource,sNormFactor,vNormFactor);

    // Calculate initial residual
    matrix_.matrixMul(sW, vW, sTmp, vTmp);

    // Approximate initial residual
    forN(nScalar,i) sTmp[i] = sSource[i] - sTmp[i];
    forN(nVector,i) vTmp[i] = vSource[i] - vTmp[i];

    forAll(sNormFactor, i)
    {
        solverPerf.sInitRes()[i] = gSumMag(sTmp[i]) / sNormFactor[i];
        solverPerf.sFinalRes()[i] = solverPerf.sInitRes()[i];
    }

    // Use vector magnitude for normalisation
    forAll(vNormFactor, i)
    {
    	solverPerf.vInitRes()[i] = cmptDivide(gSum(cmptMag(vTmp[i])),vNormFactor[i]);
        solverPerf.vFinalRes()[i] = solverPerf.vInitRes()[i];
    }

    autoPtr<coupledMatrix::smoother> smootherPtr = coupledMatrix::smoother::New
    (
        matrix_,
        controlDict_
    );

    do                                            // Outer loop
    {
        smootherPtr->smooth
        (
            sW,
			vW,
			sSource,
			vSource,
            nSweeps_
        );

        // Re-calculate the residual
		this->matrix_.matrixMul(sW, vW, sTmp, vTmp);

        forN(nScalar,j)
        {
            forAll(sTmp[j], iCell)
            {
                sTmp[j][iCell] = sSource[j][iCell] - sTmp[j][iCell];
            }
        }
        forN(nVector,j)
        {
            forAll(vTmp[j], iCell)
            {
                vTmp[j][iCell] = vSource[j][iCell] - vTmp[j][iCell];
            }
        }

        forN(nScalar,i)
        {
        	solverPerf.sFinalRes()[i] = gSumMag(sTmp[i]) / sNormFactor[i];
        }

        forN(nVector,i)
        {
            vector res = cmptDivide(gSum(cmptMag(vTmp[i])),vNormFactor[i]);

            // Zero the residual in non-solved directions
            vector::labelType validComponents = this->matrix_.mesh().solutionD(); //-1 for empty directions

            forAll(validComponents, cmpt)
            {
                if (validComponents[cmpt] == -1)
                {
                    res[cmpt] = 0.0;
                }
            }

            solverPerf.vFinalRes()[i] = res;
        }

        solverPerf.nIterations() += nSweeps_;

    } while (!this->stop(solverPerf));

    return solverPerf;

}

residualsIO smoothSolverCoupled::solveDelta
(
    PtrList<volScalarField>& sW, PtrList<volVectorField>& vW,
    const PtrList<scalarField>& sSource, const PtrList<vectorField>& vSource
) const
{
    const int nScalar = this->matrix().nScal();
    const int nVector = this->matrix().nVect();

	residualsIO solverPerf(nScalar,nVector);

    // Allocate variables to hold solved increment
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
                    "d" + sW[i].name(),
					this->matrix_.mesh().time().timeName(),
					this->matrix_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                sW[i],
                zeroGradientFvPatchScalarField::typeName
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
                    "d" + vW[i].name(),
					this->matrix_.mesh().time().timeName(),
					this->matrix_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                vW[i],
                zeroGradientFvPatchVectorField::typeName
            )
        );
    }

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

    scalarList sNormFactor(nScalar);
    scalarList vNormFactor(nVector);

    forN(nScalar,i)
    {
        scalar avg = gAverage(sW[i]);
        dsW[i].primitiveFieldRef() -= avg;
        dsW[i].boundaryFieldRef() -= avg;
    }

    forN(nVector,i)
    {
        vector avg = gAverage(vW[i]);
        dvW[i].primitiveFieldRef()-= avg;
        dvW[i].boundaryFieldRef() -= avg;
    }

    // Matrix multiplication
    this->matrix_.matrixMul(dsW, dvW, sTmp, vTmp);

    forN(nScalar,i)
    {
        dsW[i].primitiveFieldRef() = 0.0;
    }
    forN(nVector,i)
    {
        dvW[i].primitiveFieldRef() = vector::zero;
    }

    forAll(sNormFactor, i)
    {
        sNormFactor[i] = gSum(mag(sTmp[i]) + mag(sSource[i])) + VSMALL;
        solverPerf.sInitRes()[i] = gSumMag(sSource[i]) / sNormFactor[i];
        solverPerf.sFinalRes()[i] = solverPerf.sInitRes()[i];
    }

    // Use vector magnitude for normalisation
    forAll(vNormFactor, i)
    {
        vNormFactor[i] = gSum(mag(vTmp[i]) + mag(vSource[i])) + VSMALL;
        solverPerf.vInitRes()[i] = gSumCmptMag(vSource[i])/vNormFactor[i];
        solverPerf.vFinalRes()[i] = solverPerf.vInitRes()[i];
    }

    // Approximate initial residual
    forN(nScalar,i) sTmp[i] = sSource[i];
    forN(nVector,i) vTmp[i] = vSource[i];

    autoPtr<coupledMatrix::smoother> smootherPtr = coupledMatrix::smoother::New
    (
        matrix_,
        controlDict_
    );

    do                                            // Outer loop
    {
        smootherPtr->smooth
        (
            dsW,
			dvW,
			sSource,
			vSource,
            nSweeps_
        );

        // Re-calculate the residual
		this->matrix_.matrixMul(dsW, dvW, sTmp, vTmp);

        forN(nScalar,j)
        {
            forAll(sTmp[j], iCell)
            {
                sTmp[j][iCell] = sSource[j][iCell] - sTmp[j][iCell];
            }
        }
        forN(nVector,j)
        {
            forAll(vTmp[j], iCell)
            {
                vTmp[j][iCell] = vSource[j][iCell] - vTmp[j][iCell];
            }
        }

        forN(nScalar,i)
        {
        	solverPerf.sFinalRes()[i] = gSumMag(sTmp[i]) / sNormFactor[i];
        }

        forN(nVector,i)
        {
        	vector res = gSumCmptMag(vTmp[i])/vNormFactor[i];

            // Zero the residual in non-solved directions
            vector::labelType validComponents = this->matrix_.mesh().solutionD(); //-1 for empty directions

            forAll(validComponents, cmpt)
            {
                if (validComponents[cmpt] == -1)
                {
                    res[cmpt] = 0.0;
                }
            }

            solverPerf.vFinalRes()[i] = res;
        }

        solverPerf.nIterations() += nSweeps_;

    } while (!this->stop(solverPerf));

    forN(nScalar,i)
	{
		sW[i].primitiveFieldRef() += dsW[i].primitiveFieldRef();
	}
	forN(nVector,i)
	{
		vW[i].primitiveFieldRef() += dvW[i].primitiveFieldRef();
	}

    return solverPerf;

}


residualsIO smoothSolverCoupled::solveDelta
(
    PtrList<volScalarField>& sW, PtrList<volVectorField>& vW,
    const PtrList<scalarField>& sSource, const PtrList<vectorField>& vSource,
	PtrList<volScalarField>& dsW, PtrList<volVectorField>& dvW
) const
{
    const int nScalar = this->matrix().nScal();
    const int nVector = this->matrix().nVect();

	residualsIO solverPerf(nScalar,nVector);

    forN(nScalar,i)
    {
    	dsW[i].primitiveFieldRef() = sW[i].primitiveField();
    	dsW[i].boundaryFieldRef() = sW[i].boundaryField();
    }
    forN(nVector,i)
    {
    	dvW[i].primitiveFieldRef() = vW[i].primitiveField();
    	dvW[i].boundaryFieldRef() = vW[i].boundaryField();
    }

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

    scalarList sNormFactor(nScalar);
    scalarList vNormFactor(nVector);

    forN(nScalar,i)
    {
        scalar avg = gAverage(sW[i]);
        dsW[i].primitiveFieldRef() -= avg;
        dsW[i].boundaryFieldRef() -= avg;
    }
    forN(nVector,i)
    {
        vector avg = gAverage(vW[i]);
        dvW[i].primitiveFieldRef()-= avg;
        dvW[i].boundaryFieldRef() -= avg;
    }

    // Matrix multiplication
    this->matrix_.matrixMul(dsW, dvW, sTmp, vTmp);

    forN(nScalar,i)
    {
        dsW[i].primitiveFieldRef() = 0.0;
    }
    forN(nVector,i)
    {
        dvW[i].primitiveFieldRef() = vector::zero;
    }

    forAll(sNormFactor, i)
    {
        sNormFactor[i] = gSum(mag(sTmp[i]) + mag(sSource[i])) + VSMALL;
        solverPerf.sInitRes()[i] = gSumMag(sSource[i]) / sNormFactor[i];
        solverPerf.sFinalRes()[i] = solverPerf.sInitRes()[i];
    }

    // Use vector magnitude for normalisation
    forAll(vNormFactor, i)
    {
        vNormFactor[i] = gSum(mag(vTmp[i]) + mag(vSource[i])) + VSMALL;
        solverPerf.vInitRes()[i] = gSumCmptMag(vSource[i])/vNormFactor[i];
        solverPerf.vFinalRes()[i] = solverPerf.vInitRes()[i];
    }

    // Approximate initial residual
    forN(nScalar,i) sTmp[i] = sSource[i];
    forN(nVector,i) vTmp[i] = vSource[i];

    autoPtr<coupledMatrix::smoother> smootherPtr = coupledMatrix::smoother::New
    (
        matrix_,
        controlDict_
    );

    do                                            // Outer loop
    {
        smootherPtr->smooth
        (
            dsW,
			dvW,
			sSource,
			vSource,
            nSweeps_
        );

        // Re-calculate the residual
		this->matrix_.matrixMul(dsW, dvW, sTmp, vTmp);

        forN(nScalar,j)
        {
            forAll(sTmp[j], iCell)
            {
                sTmp[j][iCell] = sSource[j][iCell] - sTmp[j][iCell];
            }
        }
        forN(nVector,j)
        {
            forAll(vTmp[j], iCell)
            {
                vTmp[j][iCell] = vSource[j][iCell] - vTmp[j][iCell];
            }
        }

        forN(nScalar,i)
        {
        	solverPerf.sFinalRes()[i] = gSumMag(sTmp[i]) / sNormFactor[i];
        }

        forN(nVector,i)
        {
        	vector res = gSumCmptMag(vTmp[i])/vNormFactor[i];

            // Zero the residual in non-solved directions
            vector::labelType validComponents = this->matrix_.mesh().solutionD(); //-1 for empty directions

            forAll(validComponents, cmpt)
            {
                if (validComponents[cmpt] == -1)
                {
                    res[cmpt] = 0.0;
                }
            }

            solverPerf.vFinalRes()[i] = res;
        }

        solverPerf.nIterations() += nSweeps_;

    } while (!this->stop(solverPerf));

    return solverPerf;

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
