/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

    Hrvoje Jasak, Wikki Ltd. All rights reserved
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


#include "gmres.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

	defineTypeNameAndDebug(gmres, 0);

	coupledMatrix::solver::adddictionaryConstructorToTable<gmres>
    	addgmresDictionaryConstructorToTable_;

// * * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * //

inline void gmres::givensRotation
(
    const scalar& h,
    const scalar& beta,
    scalar& c,
    scalar& s
) const
{
    if (beta == 0)
    {
        c = 1;
        s = 0;
    }
    else if (mag(beta) > mag(h))
    {
        scalar tau = -h/beta;
        s = 1.0/Foam::sqrt(1.0 + sqr(tau));
        c = s*tau;
    }
    else
    {
        scalar tau = -beta/h;
        c = 1.0/Foam::sqrt(1.0 + sqr(tau));
        s = c*tau;
    }
}


// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

gmres::gmres
(
    const dictionary& dict,
    const coupledMatrix& matrix
)
:
	coupledMatrix::solver(typeName, dict, matrix),
	nDirs_(readLabel(dict.lookup("nDirections")))
{}


residualsIO gmres::solve
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

    // Create the Hessenberg matrix
    scalarSquareMatrix H(nDirs_, Zero);	        // Initialise H [m x m]

    // Create y and b for Hessenberg matrix
    scalarField yh(nDirs_, 0);
    scalarField bh(nDirs_ + 1, 0);

    // Givens rotation vectors
    scalarField c(nDirs_, 0);
    scalarField s(nDirs_, 0);

    // Allocate Krylov vectors
    List<PtrList<volScalarField > > sVPtr(nDirs_);
    List<PtrList<volVectorField > > vVPtr(nDirs_);

    for(label i = 0; i < nDirs_; i++)
    {
        sVPtr[i].setSize(nScalar);

        forN(nScalar,j)
        {
            sVPtr[i].set
            (
                j,
                new volScalarField
                (
                    IOobject
                    (
                        "krylovVectorScalars" + Foam::name(i) + "-" + Foam::name(j),
						matrix_.mesh().time().timeName(),
						matrix_.mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
					matrix_.mesh(),
                    dimless,
                    zeroGradientFvPatchScalarField::typeName
                )
            );
        }

        vVPtr[i].setSize(nVector);

        forN(nVector,j)
        {
            vVPtr[i].set
            (
                j,
                new volVectorField
                (
                    IOobject
                    (
                        "krylovVectorVectors" + Foam::name(i) + "-" + Foam::name(j),
						matrix_.mesh().time().timeName(),
						matrix_.mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    matrix_.mesh(),
                    dimless,
                    zeroGradientFvPatchVectorField::typeName
                )
            );
        }
    }

    // --- Select and construct the preconditioner
    autoPtr<coupledMatrix::preconditioner> preconPtr =
        coupledMatrix::preconditioner::New
        (
            *this,
            controlDict_
        );

    do                                            // Outer loop
    {
        // Execute preconditioning
        preconPtr->precondition(sTmp, vTmp);  // P^-1 v_o  &  P^-1 A \Delta U_0

        // Calculate beta
        scalar beta = 0.0;
        forN(nScalar,i)
        {
            beta += gSumSqr(sTmp[i]);
        }
        forN(nVector,i)
        {
            beta += gSum(magSqr(vTmp[i]));
        }

        beta = Foam::sqrt(beta);

        // Set initial rhs and bh[0] = beta
        bh = 0;
        bh[0] = beta;

        for (label i = 0; i < nDirs_; i++)				// Search directions
        {
            // Set search direction - delay scaling vector to allow parallel comms overlap
            PtrList<volScalarField>& sV = sVPtr[i];
            PtrList<volVectorField>& vV = vVPtr[i];

            forN(nScalar,j)
            {
                sV[j].primitiveFieldRef() = sTmp[j]/beta;
            }
            forN(nVector,j)
            {
                vV[j].primitiveFieldRef() = vTmp[j]/beta;
            }

            // Matrix vector product
            this->matrix_.matrixMul(sV, vV, sTmp, vTmp);    // y_j = A v_j

            // Execute preconditioning
            preconPtr->precondition(sTmp, vTmp);			// w_j = P^-1 y_j

            for(label j = 0; j <= i; j++)
            {
                beta = 0.0;

                forN(nScalar,k)
                {
                    beta += gSumProd(sTmp[k], sVPtr[j][k]);
                }

                forN(nVector,k)
                {
                    beta += gSumProd(vTmp[k], vVPtr[j][k]);
                }

            	H[j][i] = beta;

                forN(nScalar,k)                          // w_j = w_j - h_ij v_i
                {
                    sTmp[k] -= H[j][i]*sVPtr[j][k].primitiveFieldRef();
                }
                forN(nVector,k)
                {
                    vTmp[k] -= H[j][i]*vVPtr[j][k].primitiveFieldRef();
                }
            }

            beta = 0.0;

            forN(nScalar,j)
            {
                beta += gSumSqr(sTmp[j]);
            }

            forN(nVector,j)
            {
                beta += gSum(magSqr(vTmp[j]));
            }

            beta = Foam::sqrt(beta);

            // Apply previous Givens rotations to new column of H.
            for (label j = 0; j < i; j++)
            {
                const scalar Hji = H[j][i];				// Givens rotation similar to Saad
                H[j][i] = c[j]*Hji - s[j]*H[j+1][i];
                H[j+1][i] = s[j]*Hji + c[j]*H[j+1][i];
            }

            // Apply Givens rotation to current row
            givensRotation(H[i][i], beta, c[i], s[i]);

            const scalar bhi = bh[i];
            bh[i] = c[i]*bhi - s[i]*bh[i+1];
            bh[i + 1] = s[i]*bhi + c[i]*bh[i+1];
            H[i][i] = c[i]*H[i][i] - s[i]*beta;

        }

        // Back substitute to solve Hy = b
        for (label i = nDirs_ - 1; i >= 0; i--)
        {
            scalar sum = bh[i];

            for (label j = i + 1; j < nDirs_; j++)
            {
                sum -= H[i][j]*yh[j];
            }

            yh[i] = sum/stabilise(H[i][i], VSMALL);  // In case of zero initial residual
        }

        // Update solution
        for (label i = 0; i < nDirs_; i++)
        {
            const PtrList<volScalarField>& sVi = sVPtr[i];
            const PtrList<volVectorField>& vVi = vVPtr[i];

            const scalar& yi = yh[i];

            forN(nScalar,j)
            {
            	sW[j].primitiveFieldRef() += yi*sVi[j].primitiveField();
            }
            forN(nVector,j)
            {
                vW[j].primitiveFieldRef() += yi*vVi[j].primitiveField();
            }
        }

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

        solverPerf.nIterations()++;

    } while (!this->stop(solverPerf));

    return solverPerf;

}

residualsIO gmres::solveDelta
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

    // Create the Hessenberg matrix
    scalarSquareMatrix H(nDirs_, Zero);	        // Initialise H [m x m]

    // Create y and b for Hessenberg matrix
    scalarField yh(nDirs_, 0);
    scalarField bh(nDirs_ + 1, 0);

    // Givens rotation vectors
    scalarField c(nDirs_, 0);
    scalarField s(nDirs_, 0);

    // Allocate Krylov vectors
     List<PtrList<volScalarField > > sVPtr(nDirs_);
     List<PtrList<volVectorField > > vVPtr(nDirs_);

     for(label i = 0; i < nDirs_; i++)
     {
         sVPtr[i].setSize(nScalar);

         forN(nScalar,j)
         {
             sVPtr[i].set
             (
                 j,
                 new volScalarField
                 (
                     IOobject
                     (
						 "krylovVectorScalars" + Foam::name(i) + "-" + Foam::name(j),
 						 matrix_.mesh().time().timeName(),
 						 matrix_.mesh(),
                         IOobject::NO_READ,
                         IOobject::NO_WRITE
                     ),
 					matrix_.mesh(),
                     dimless,
                     zeroGradientFvPatchScalarField::typeName
                 )
             );
         }

         vVPtr[i].setSize(nVector);

         forN(nVector,j)
         {
             vVPtr[i].set
             (
                 j,
                 new volVectorField
                 (
                     IOobject
                     (
						 "krylovVectorVectors" + Foam::name(i) + "-" + Foam::name(j),
 						matrix_.mesh().time().timeName(),
 						matrix_.mesh(),
                         IOobject::NO_READ,
                         IOobject::NO_WRITE
                     ),
                     matrix_.mesh(),
                     dimless,
                     zeroGradientFvPatchVectorField::typeName
                 )
             );
         }
     }

    // --- Select and construct the preconditioner
    autoPtr<coupledMatrix::preconditioner> preconPtr =
        coupledMatrix::preconditioner::New
        (
            *this,
            controlDict_
        );

    do                                            // Outer loop
    {
        // Execute preconditioning
        preconPtr->precondition(sTmp, vTmp);  // P^-1 v_o  &  P^-1 A \Delta U_0

        // Calculate beta
        scalar beta = 0.0;
        forN(nScalar,i)
        {
            beta += gSumSqr(sTmp[i]);
        }
        forN(nVector,i)
        {
            beta += gSum(magSqr(vTmp[i]));
        }

        beta = Foam::sqrt(beta);

        // Set initial rhs and bh[0] = beta
        bh = 0;
        bh[0] = beta;

        for (label i = 0; i < nDirs_; i++)				// Search directions
        {
            // Set search direction - delay scaling vector to allow parallel comms overlap
            PtrList<volScalarField>& sV = sVPtr[i];
            PtrList<volVectorField>& vV = vVPtr[i];

            forN(nScalar,j)
            {
                sV[j].primitiveFieldRef() = sTmp[j]/beta;
            }
            forN(nVector,j)
            {
                vV[j].primitiveFieldRef() = vTmp[j]/beta;
            }

            // Matrix vector product
            this->matrix_.matrixMul(sV, vV, sTmp, vTmp);    // y_j = A v_j

            // Execute preconditioning
            preconPtr->precondition(sTmp, vTmp);			// w_j = P^-1 y_j

            for(label j = 0; j <= i; j++)
            {
                beta = 0.0;

                forN(nScalar,k)
                {
                    beta += gSumProd(sTmp[k], sVPtr[j][k]);
                }

                forN(nVector,k)
                {
                    beta += gSumProd(vTmp[k], vVPtr[j][k]);
                }

            	H[j][i] = beta;

                forN(nScalar,k)                          // w_j = w_j - h_ij v_i
                {
                    sTmp[k] -= H[j][i]*sVPtr[j][k];
                }
                forN(nVector,k)
                {
                    vTmp[k] -= H[j][i]*vVPtr[j][k];
                }
            }

            beta = 0.0;

            forN(nScalar,j)
            {
                beta += gSumSqr(sTmp[j]);
            }

            forN(nVector,j)
            {
                beta += gSum(magSqr(vTmp[j]));
            }

            beta = Foam::sqrt(beta);

            // Apply previous Givens rotations to new column of H.
            for (label j = 0; j < i; j++)
            {
                const scalar Hji = H[j][i];				// Givens rotation similar to Saad
                H[j][i] = c[j]*Hji - s[j]*H[j+1][i];
                H[j+1][i] = s[j]*Hji + c[j]*H[j+1][i];
            }

            // Apply Givens rotation to current row
            givensRotation(H[i][i], beta, c[i], s[i]);

            const scalar bhi = bh[i];
            bh[i] = c[i]*bhi - s[i]*bh[i+1];
            bh[i + 1] = s[i]*bhi + c[i]*bh[i+1];
            H[i][i] = c[i]*H[i][i] - s[i]*beta;

        }

        // Back substitute to solve Hy = b
        for (label i = nDirs_ - 1; i >= 0; i--)
        {
            scalar sum = bh[i];

            for (label j = i + 1; j < nDirs_; j++)
            {
                sum -= H[i][j]*yh[j];
            }

            yh[i] = sum/stabilise(H[i][i], VSMALL);  // In case of zero initial residual
        }

        // Update solution
        for (label i = 0; i < nDirs_; i++)
        {
            const PtrList<volScalarField>& sVi = sVPtr[i];
            const PtrList<volVectorField>& vVi = vVPtr[i];

            const scalar& yi = yh[i];

            forN(nScalar,j)
            {
            	dsW[j].primitiveFieldRef() += yi*sVi[j].primitiveField();
            }
            forN(nVector,j)
            {
                dvW[j].primitiveFieldRef() += yi*vVi[j].primitiveField();
            }
        }

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

        solverPerf.nIterations()++;

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


residualsIO gmres::solveDelta
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

    // Create the Hessenberg matrix
    scalarSquareMatrix H(nDirs_, Zero);	        // Initialise H [m x m]

    // Create y and b for Hessenberg matrix
    scalarField yh(nDirs_, 0);
    scalarField bh(nDirs_ + 1, 0);

    // Givens rotation vectors
    scalarField c(nDirs_, 0);
    scalarField s(nDirs_, 0);

    // Allocate Krylov vectors
     List<PtrList<volScalarField > > sVPtr(nDirs_);
     List<PtrList<volVectorField > > vVPtr(nDirs_);

     for(label i = 0; i < nDirs_; i++)
     {
         sVPtr[i].setSize(nScalar);

         forN(nScalar,j)
         {
             sVPtr[i].set
             (
                 j,
                 new volScalarField
                 (
                     IOobject
                     (
						 "krylovVectorScalars" + Foam::name(i) + "-" + Foam::name(j),
 						 sW[j].mesh().time().timeName(),
 						 sW[j].mesh(),
                         IOobject::NO_READ,
                         IOobject::NO_WRITE
                     ),
 					 sW[j].mesh(),
                     dimless,
                     zeroGradientFvPatchScalarField::typeName
                 )
             );
         }

         vVPtr[i].setSize(nVector);

         forN(nVector,j)
         {
             vVPtr[i].set
             (
                 j,
                 new volVectorField
                 (
                     IOobject
                     (
						 "krylovVectorVectors" + Foam::name(i) + "-" + Foam::name(j),
 						 vW[j].mesh().time().timeName(),
 						 vW[j].mesh(),
                         IOobject::NO_READ,
                         IOobject::NO_WRITE
                     ),
                     vW[j].mesh(),
                     dimless,
                     zeroGradientFvPatchVectorField::typeName
                 )
             );
         }
     }

    // --- Select and construct the preconditioner
    autoPtr<coupledMatrix::preconditioner> preconPtr =
        coupledMatrix::preconditioner::New
        (
            *this,
            controlDict_
        );

    do                                            // Outer loop
    {
        // Execute preconditioning
        preconPtr->precondition(sTmp, vTmp);  // P^-1 v_o  &  P^-1 A \Delta U_0

        // Calculate beta
        scalar beta = 0.0;
        forN(nScalar,i)
        {
            beta += gSumSqr(sTmp[i]);
        }
        forN(nVector,i)
        {
            beta += gSum(magSqr(vTmp[i]));
        }

        beta = Foam::sqrt(beta);

        // Set initial rhs and bh[0] = beta
        bh = 0;
        bh[0] = beta;

        for (label i = 0; i < nDirs_; i++)				// Search directions
        {
            // Set search direction - delay scaling vector to allow parallel comms overlap
            PtrList<volScalarField>& sV = sVPtr[i];
            PtrList<volVectorField>& vV = vVPtr[i];

            forN(nScalar,j)
            {
                sV[j].primitiveFieldRef() = sTmp[j]/beta;
            }
            forN(nVector,j)
            {
                vV[j].primitiveFieldRef() = vTmp[j]/beta;
            }

            // Matrix vector product
            this->matrix_.matrixMul(sV, vV, sTmp, vTmp);    // y_j = A v_j

            // Execute preconditioning
            preconPtr->precondition(sTmp, vTmp);			// w_j = P^-1 y_j

            for(label j = 0; j <= i; j++)
            {
                beta = 0.0;

                forN(nScalar,k)
                {
                    beta += gSumProd(sTmp[k], sVPtr[j][k]);
                }

                forN(nVector,k)
                {
                    beta += gSumProd(vTmp[k], vVPtr[j][k]);
                }

            	H[j][i] = beta;

                forN(nScalar,k)                          // w_j = w_j - h_ij v_i
                {
                    sTmp[k] -= H[j][i]*sVPtr[j][k];
                }
                forN(nVector,k)
                {
                    vTmp[k] -= H[j][i]*vVPtr[j][k];
                }
            }

            beta = 0.0;

            forN(nScalar,j)
            {
                beta += gSumSqr(sTmp[j]);
            }

            forN(nVector,j)
            {
                beta += gSum(magSqr(vTmp[j]));
            }

            beta = Foam::sqrt(beta);

            // Apply previous Givens rotations to new column of H.
            for (label j = 0; j < i; j++)
            {
                const scalar Hji = H[j][i];				// Givens rotation similar to Saad
                H[j][i] = c[j]*Hji - s[j]*H[j+1][i];
                H[j+1][i] = s[j]*Hji + c[j]*H[j+1][i];
            }

            // Apply Givens rotation to current row
            givensRotation(H[i][i], beta, c[i], s[i]);

            const scalar bhi = bh[i];
            bh[i] = c[i]*bhi - s[i]*bh[i+1];
            bh[i + 1] = s[i]*bhi + c[i]*bh[i+1];
            H[i][i] = c[i]*H[i][i] - s[i]*beta;
        }

        // Back substitute to solve Hy = b
        for (label i = nDirs_ - 1; i >= 0; i--)
        {
            scalar sum = bh[i];

            for (label j = i + 1; j < nDirs_; j++)
            {
                sum -= H[i][j]*yh[j];
            }

            yh[i] = sum/stabilise(H[i][i], VSMALL);  // In case of zero initial residual
        }

        // Update solution
        for (label i = 0; i < nDirs_; i++)
        {
            const PtrList<volScalarField>& sVi = sVPtr[i];
            const PtrList<volVectorField>& vVi = vVPtr[i];

            const scalar& yi = yh[i];

            forN(nScalar,j)
            {
            	dsW[j].primitiveFieldRef() += yi*sVi[j].primitiveField();
            }
            forN(nVector,j)
            {
                dvW[j].primitiveFieldRef() += yi*vVi[j].primitiveField();
            }
        }

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

        solverPerf.nIterations()++;

    } while (!this->stop(solverPerf));

    return solverPerf;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
