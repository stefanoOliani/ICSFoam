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

#include "HBZone.H"

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"
#include "faceSet.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "EigenMatrix.H"
#include "LUscalarMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HBZone, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::HBZone::setOmegaList()
{
	scalarList freqList = coeffs_.get<scalarList>("frequenciesList");
	labelList harmList = coeffs_.get<labelList>("harmonicsList");

	if (freqList.size() != harmList.size())
	{
		FatalErrorInFunction
			<< "harmonicsList size must match frequenciesList size "
			<< exit(FatalError);
	}

	// The first frequency is always zero (steady term)
	omegaList_.append(0.0);

	scalarList posFrequencies;
	scalarList negFrequencies;

	forAll (freqList,i)
	{
		scalar freqi = freqList[i];

		//- Always work with positive frequencies
		if (freqi < 0.0)
		{
			freqi = -freqi;
		}

		label harmi = harmList[i];

		for (int k = 1; k <= harmi; k++)
		{
			posFrequencies.append(k*freqi);
			negFrequencies.append(-k*freqi);
		}
	}

	omegaList_.append(posFrequencies);

	forAllReverse (negFrequencies, i)
	{
		omegaList_.append(negFrequencies[i]);
	}

	Info<<"Omega list "<<omegaList_<<endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HBZone::HBZone
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    mesh_(mesh),
    name_(name),
    coeffs_(dict),
    active_(coeffs_.getOrDefault("active", true)),
    cellZoneName_(cellZoneName),
    cellZoneID_(),
	omegaList_(),
	cylCoords_(dict.get<bool>("cylCoords")),
	rotationAxis_(Zero),
	rotationCentre_(Zero),
	EInv_(),
	E_(),
	D_()
{
    if (cellZoneName_ == word::null)
    {
        // If the source term of this zone has to be
        // applied to the whole mesh
        if (coeffs_.getOrDefault("allMesh", false))
        {
        	cellZoneID_ = -2;
        }
        else
        {
            coeffs_.readEntry("cellZone", cellZoneName_);
        }
    }

    if (!active_)
    {
        cellZoneID_ = -1;
    }
    else if (cellZoneID_ == -2)
	{
    	//Do nothing
	}
    else
    {
        cellZoneID_ = mesh_.cellZones().findZoneID(cellZoneName_);

        bool cellZoneFound = (cellZoneID_ != -1);

        reduce(cellZoneFound, orOp<bool>());

        if (!cellZoneFound)
        {
            FatalErrorInFunction
                << "cannot find HB cellZone " << cellZoneName_
                << exit(FatalError);
        }
    }

    this->setOmegaList();

    const dictionary& solDict = mesh_.solutionDict();

	const label nInstants =
			solDict.subDict("harmonicBalance").getOrDefault<label>("instantsNumber",3);

	bool oversampling =
			solDict.subDict("harmonicBalance").getOrDefault<bool>("oversampling",false);

	Info<<"For zone "<<name<<" the cell zone ID is "<<cellZoneID_<<endl;

	if (nInstants != omegaList_.size())
	{
		if(!oversampling)
		{
			FatalErrorInFunction
				<< "specified number of instants is not correct "
				<< exit(FatalError);
		}

		if(nInstants < omegaList_.size())
		{
			FatalErrorInFunction
				<< "specified number of instants is lower than frequency vector dimension "
				<< exit(FatalError);
		}

	}

	if (cylCoords_)
	{
		coeffs_.readEntry("rotationAxis", rotationAxis_);
		coeffs_.readEntry("rotationCentre", rotationCentre_);
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::HBZone::calcConditionN(const scalarList& snapshots)
{
	if (snapshots.size() < omegaList_.size())
	{
		FatalErrorInFunction
			<< "snapshots size must be equal or greater to omegaList size "
			<< exit(FatalError);
	}

	RectangularMatrix<complex> E_1(snapshots.size(),omegaList_.size());
	complex e(0,1);

	forAll (snapshots , n)
	{
		forAll (omegaList_, k)
		{
			e.Re() = Foam::cos(omegaList_[k]*snapshots[n]);
			e.Im() = Foam::sin(omegaList_[k]*snapshots[n]);
			E_1[n][k] = e;
		}
	}

	//- Structure to compute L2 norm of the matrices

	RectangularMatrix<complex> E_1TrConj = E_1.T();

	SquareMatrix<complex> E_1Mult(omegaList_.size(), Zero);
	SquareMatrix<scalar> E_1MultRe(omegaList_.size(), Zero);


	for (label l=0; l<omegaList_.size(); l++)
	{
		for (label m=0; m<omegaList_.size(); m++)
		{
			for (label n=0; n<snapshots.size(); n++)
			{
				E_1Mult[l][m] += E_1TrConj[l][n]*E_1[n][m];
				E_1MultRe[l][m] = E_1Mult[l][m].Re();
			}
		}
	}

	//Info<<"E_1Mult "<<E_1Mult<<endl;

	//- Now we compute the SVD of matrix E_1 resorting to
	//  the corresponding eigenvalue problem E_1E_1*
	//  (this is necessary because the SVD class of OF only
	//  takes real matrices as input)

	EigenMatrix<scalar> E_1SVD(E_1MultRe);

	//- Singular values are the square root of the
	//  eigenvalues of E_1E_1*. Moreover, we are sure that
	//  the eigenvalues are non-negative real numbers and
	//  therefore  it suffices to use the EValsRe()
	//  function of theEigenMatrix class.

	List<scalar> singularValues = E_1SVD.EValsRe();

	if (min(singularValues) < SMALL)
	{
		return GREAT;
	}

	scalar maxSV = Foam::sqrt(max(singularValues));
	scalar minSV = Foam::sqrt(min(singularValues));

	//- The L2 norm condition number of a
	//  rectangular matrix is defined as the ratio
	//  between the maximum and minimum singular value

	return maxSV/minSV;

}

void Foam::HBZone::updateHBOperators(const scalarList& selSnapshots)
{
	const label nO = selSnapshots.size(); //number of snapshots
	const label nT = omegaList_.size(); //number of frequencies

	RectangularMatrix<complex> E_1(nO, nT);
	complex e(0,1);

	forAll (selSnapshots , n)
	{
		forAll (omegaList_, k)
		{
			e.Re() = Foam::cos(omegaList_[k]*selSnapshots[n]);
			e.Im() = Foam::sin(omegaList_[k]*selSnapshots[n]);
			E_1[n][k] = e;
		}
	}

	EInv_ = E_1;

	Info<<"EInv "<<EInv_<<endl;

	E_ = Foam::pinv(E_1);

	Info<<"E "<<E_<<endl;

	// Decompose the transform matrices into
	// real and imaginary part. D can be
	// computed as -(EInvIm*A*ERe + EInvRe*A*EIm)

	SquareMatrix<scalar> A(nT, 0.0);

	RectangularMatrix<scalar> ERe(nT, nO, 0.0);
	RectangularMatrix<scalar> EIm(nT, nO, 0.0);
	RectangularMatrix<scalar> EInvRe(nO, nT, 0.0);
	RectangularMatrix<scalar> EInvIm(nO, nT, 0.0);

	for (label l=0; l<nT; l++)
	{
		A[l][l] = omegaList_[l];

		for (label m=0; m<nO; m++)
		{
			ERe[l][m] = E_[l][m].Re();
			EIm[l][m] = E_[l][m].Im();
			EInvRe[m][l] = EInv_[m][l].Re();
			EInvIm[m][l] = EInv_[m][l].Im();
		}
	}

	// Matrices defined just to easy the
	// computation of D in several steps

	RectangularMatrix<scalar> temp_0(nO, nT, 0.0);
	RectangularMatrix<scalar> temp_1(nO, nT, 0.0);

	SquareMatrix<scalar> d_0(nO, 0.0);
	SquareMatrix<scalar> d_1(nO, 0.0);

	for (label l=0; l<nO; l++)
	{
		for (label m=0; m<nT; m++)
		{
			for (label n=0; n<nT; n++)
			{
				temp_0[l][m] += EInvIm[l][n]*A[n][m];
				temp_1[l][m] += EInvRe[l][n]*A[n][m];
			}
		}
	}

	for (label l=0; l<nO; l++)
	{
		for (label m=0; m<nO; m++)
		{
			for (label n=0; n<nT; n++)
			{
				d_0[l][m] += temp_0[l][n]*ERe[n][m];
				d_1[l][m] += temp_1[l][n]*EIm[n][m];
			}
		}
	}

	D_ = -(d_0 + d_1);

	Info<<"D matrix "<<D_<<endl;
}


void Foam::HBZone::addSource
(
	coupledMatrix& eqSystem,
	PtrList<PtrList<volScalarField>>& scalarVars,
	PtrList<PtrList<volVectorField>>& vectorVars,
	label instantNo
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const scalarField& V = mesh_.V();

    if (cellZoneID_ == -2)
    {
    	const labelList cells(mesh_.nCells());

    	label cellNo = mesh_.nCells();

    	for (Foam::label celli=0; celli<cellNo; ++celli)
    	{
    		forAll (scalarVars[instantNo], i)
			{
    			scalarField& scalSource = eqSystem.dSByS(i,i).source();

        		forAll (scalarVars, J)
        		{
        			scalSource[celli] -= V[celli]*this->D()[instantNo][J]*scalarVars[J][i][celli];
        		}
			}

    		forAll (vectorVars[instantNo], i)
			{
				vectorField& vectorSource = eqSystem.dVByV(i,i).source();

				forAll (vectorVars, J)
				{
					vectorSource[celli] -= V[celli]*this->D()[instantNo][J]*vectorVars[J][i][celli];
				}
			}
    	}
    }
    else
    {
    	const labelList& cells = mesh_.cellZones()[cellZoneID_];

    	forAll (cells, index)
    	{
			label celli = cells[index];

			forAll (scalarVars[instantNo], i)
			{
				scalarField& scalSource = eqSystem.dSByS(i,i).source();

				forAll (scalarVars, J)
				{
					scalSource[celli] -= V[celli]*this->D()[instantNo][J]*scalarVars[J][i][celli];
				}
			}

			forAll (vectorVars[instantNo], i)
			{
				vectorField& vectorSource = eqSystem.dVByV(i,i).source();

				forAll (vectorVars, J)
				{
					vectorSource[celli] -= V[celli]*this->D()[instantNo][J]*vectorVars[J][i][celli];
				}
			}
    	}
    }
}


void Foam::HBZone::addBlock
(
	blockFvMatrix<scalar,scalar>& blockMatrix,
	label rowIndex,
	label colIndex
) const
{
	scalarField& scalDiag = blockMatrix.diag();

    if (cellZoneID_ == -1)
    {
        return;
    }

    const scalarField& V = mesh_.V();

    if (cellZoneID_ == -2)
    {
    	const labelList cells(mesh_.nCells());

    	label cellNo = mesh_.nCells();

    	for (Foam::label celli=0; celli<cellNo; ++celli)
    	{
			//ADD because we are on the RHS
			scalDiag[celli] += V[celli]*this->D()[rowIndex][colIndex];
    	}
    }
    else
    {
    	const labelList& cells = mesh_.cellZones()[cellZoneID_];

    	forAll (cells, i)
    	{
			label celli = cells[i];

			//ADD because we are on the RHS
			scalDiag[celli] += V[celli]*this->D()[rowIndex][colIndex];
    	}
    }
}


void Foam::HBZone::addBlock
(
	blockFvMatrix<vector,tensor>& blockMatrix,
	label rowIndex,
	label colIndex
) const
{
	tensorField& tensDiag = blockMatrix.diag();

    if (cellZoneID_ == -1)
    {
        return;
    }

    const scalarField& V = mesh_.V();

    if (cellZoneID_ == -2)
    {
    	const labelList cells(mesh_.nCells());

    	label cellNo = mesh_.nCells();

    	for (Foam::label celli=0; celli<cellNo; ++celli)
    	{
			//ADD because we are on the RHS
			tensDiag[celli] += V[celli]*this->D()[rowIndex][colIndex]*tensor::I;
    	}
    }
    else
    {
    	const labelList& cells = mesh_.cellZones()[cellZoneID_];

    	forAll (cells, i)
    	{
			label celli = cells[i];

			//ADD because we are on the RHS
			tensDiag[celli] += V[celli]*this->D()[rowIndex][colIndex]*tensor::I;
    	}
    }
}


void Foam::HBZone::addSource
(
	PtrList<vectorField>& source,
	PtrList<PtrList<volVectorField>>& vars,
	label varNum,
	bool tryCylCoords
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    if (cylCoords_)
    {
    	dimensionedScalar smallRadius("smallRadius", dimLength, SMALL);
		const scalarField& V = mesh_.V();

		if (cellZoneID_ == -2)
		{
			label cellNo = mesh_.nCells();

			vector axisHat = rotationAxis_/mag(rotationAxis_);

			forAll(source, J)
			{
				source[J] = Zero;
			}

			for (Foam::label celli=0; celli<cellNo; ++celli)
			{
				forAll (source, J)
				{
					vector sourceCyl(Zero);

					forAll (source, K)
					{
						const fvMesh& meshJ = vars[K][varNum].mesh();

						const vector& cellCentre = meshJ.C()[celli];

						 // Radius vector in plane of rotation
						vector r(cellCentre - rotationCentre_);
						r -= (axisHat & r)*axisHat;
						const scalar magr(mag(r));
						const vector rHat(r/magr);

						scalar Ur = vars[K][varNum][celli] & rHat;
						scalar Uu = vars[K][varNum][celli] & (axisHat ^ rHat);
						scalar Uz = vars[K][varNum][celli] & axisHat;

						vector UCyl(Ur, Uu, Uz);

						sourceCyl += V[celli]*this->D()[J][K]*UCyl;
					}

					const vector& cellCentreAct = vars[J][varNum].mesh().C()[celli];

					// Radius vector in plane of rotation
					vector r(cellCentreAct - rotationCentre_);
					r -= (axisHat & r)*axisHat;
					const scalar magr(mag(r));
					const vector rHat(r/magr);

					vector sourceCart = sourceCyl.x()*rHat
									  + sourceCyl.y()*(axisHat^rHat)
									  + sourceCyl.z()*axisHat;

				   source[J][celli] -= sourceCart;
				}
			}
		}
		else
		{
			const labelList& cells = mesh_.cellZones()[cellZoneID_];

			vector axisHat = rotationAxis_/mag(rotationAxis_);

			forAll (cells, i)
			{
				label celli = cells[i];

				forAll (source, J)
				{
					source[J][celli] = Zero;

					vector sourceCyl(Zero);

					forAll (source, K)
					{
						const fvMesh& meshJ = vars[K][varNum].mesh();

						const vector& cellCentre = meshJ.C()[celli];

						 // Radius vector in plane of rotation
						vector r(cellCentre - rotationCentre_);
						r -= (axisHat & r)*axisHat;
						const scalar magr(mag(r));
						const vector rHat(r/magr);

						scalar Ur = vars[K][varNum][celli] & rHat;
						scalar Uu = vars[K][varNum][celli] & (axisHat ^ rHat);
						scalar Uz = vars[K][varNum][celli] & axisHat;

						vector UCyl(Ur, Uu, Uz);

						sourceCyl += V[celli]*this->D()[J][K]*UCyl;
					}

					const vector& cellCentreAct = vars[J][varNum].mesh().C()[celli];

					// Radius vector in plane of rotation
					vector r(cellCentreAct - rotationCentre_);
					r -= (axisHat & r)*axisHat;
					const scalar magr(mag(r));
					const vector rHat(r/magr);

					vector sourceCart = sourceCyl.x()*rHat
									  + sourceCyl.y()*(axisHat^rHat)
									  + sourceCyl.z()*axisHat;

				   source[J][celli] -= sourceCart;
				}
			}
		}
    }
    else
    {
    	addSource(source, vars, varNum);
    }
}


void Foam::HBZone::addSource
(
	fvVectorMatrix& eqn,
	PtrList<volVectorField>& fieldPtr,
	label instantNo,
	bool  cylCoords
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const scalarField& V = mesh_.V();
    vectorField& Source = eqn.source();

    dimensionedScalar smallRadius("smallRadius", dimLength, SMALL);

    if (cellZoneID_ == -2)
    {
    	const labelList cells(mesh_.nCells());

    	label cellNo = mesh_.nCells();

    	for (Foam::label celli=0; celli<cellNo; ++celli)
    	{
    		vector sourceCyl(Zero);

    		forAll (fieldPtr, J)
    		{
    			const fvMesh& meshJ = fieldPtr[J].mesh();
    		    scalar coordsX = meshJ.C()[celli].x();
    		    scalar coordsY = meshJ.C()[celli].y();
    		    scalar radius = std::max(sqrt(sqr(coordsX)+sqr(coordsY)),smallRadius.value());

    		    // Compute cylindrical ones from cartesian velocity components
    		   scalar Ur = (fieldPtr[J][celli].x()*coordsX+fieldPtr[J][celli].y()*coordsY)/radius;
    		   scalar Uu = (fieldPtr[J][celli].y()*coordsX-fieldPtr[J][celli].x()*coordsY)/radius;
    		   scalar Uz = fieldPtr[J][celli].z();

    		   vector UCyl(Ur, Uu, Uz);

    		   sourceCyl += V[celli]*this->D()[instantNo][J]*UCyl;
    		}

 		   vector sourceCart(0.0, 0.0, sourceCyl.z());

		   scalar actCoordsX = fieldPtr[instantNo].mesh().C()[celli].x();
		   scalar actCoordsY = fieldPtr[instantNo].mesh().C()[celli].y();

		   scalar actRadius = std::max(sqrt(sqr(actCoordsX)+sqr(actCoordsY)),smallRadius.value());

 	       sourceCart.replace(vector::X,(sourceCyl.x()*actCoordsX - sourceCyl.y()*actCoordsY)/actRadius);
 	       sourceCart.replace(vector::Y,(sourceCyl.x()*actCoordsY + sourceCyl.y()*actCoordsX)/actRadius);

 	       Source[celli] -= sourceCart;
    	}
    }
    else
    {
    	const labelList& cells = mesh_.cellZones()[cellZoneID_];

    	forAll (cells, i)
    	{
			label celli = cells[i];

			vector sourceCyl(Zero);

			forAll (fieldPtr, J)
			{
				const fvMesh& meshJ = fieldPtr[J].mesh();
				scalar coordsX = meshJ.C()[celli].x();
				scalar coordsY = meshJ.C()[celli].y();
				scalar radius = std::max(sqrt(sqr(coordsX)+sqr(coordsY)),smallRadius.value());

				// Compute cylindrical ones from cartesian velocity components
			   scalar Ur = (fieldPtr[J][celli].x()*coordsX+fieldPtr[J][celli].y()*coordsY)/radius;
			   scalar Uu = (fieldPtr[J][celli].y()*coordsX-fieldPtr[J][celli].x()*coordsY)/radius;
			   scalar Uz = fieldPtr[J][celli].z();

			   vector UCyl(Ur, Uu, Uz);

			   sourceCyl += V[celli]*this->D()[instantNo][J]*UCyl;
			}

		   vector sourceCart(0.0, 0.0, sourceCyl.z());

		   scalar actCoordsX = fieldPtr[instantNo].mesh().C()[celli].x();
		   scalar actCoordsY = fieldPtr[instantNo].mesh().C()[celli].y();
		   scalar actRadius = std::max(sqrt(sqr(actCoordsX)+sqr(actCoordsY)),smallRadius.value());

		   sourceCart.replace(vector::X,(sourceCyl.x()*actCoordsX - sourceCyl.y()*actCoordsY)/actRadius);
		   sourceCart.replace(vector::Y,(sourceCyl.x()*actCoordsY + sourceCyl.y()*actCoordsX)/actRadius);

		   Source[celli] -= sourceCart;
    	}
    }
}


void Foam::HBZone::factorizationStep
(
	PtrList<PtrList<volScalarField>>& scalarVarsIncr,
	PtrList<PtrList<volVectorField>>& vectorVarsIncr,
	PtrList<volScalarField>& rDeltaTField
) const
{
	label NT = scalarVarsIncr.size();

    if (cellZoneID_ == -1)
    {
        return;
    }

    if (cellZoneID_ == -2)
    {
    	label cellNo = mesh_.nCells();

    	for (Foam::label celli=0; celli<cellNo; ++celli)
    	{
    		// Compute the factorization matrix F

    		scalarList cellDeltaT(NT,0.0);
    		SquareMatrix<scalar> dTau(NT,0.0);

    		forAll(cellDeltaT, meshi)
    		{
    			cellDeltaT[meshi] = 1/rDeltaTField[meshi][celli];
				dTau[meshi][meshi] = cellDeltaT[meshi];
    		}

    		SquareMatrix<scalar> II(NT,0.0);
    		for (label l=0; l<NT; l++)
    		{
    			II[l][l] = 1.0;
    		}

    		SquareMatrix<scalar> f_0 = II + (dTau*this->D());
    		LUscalarMatrix f(f_0);
    		scalarSquareMatrix F(NT,0.0);
    		f.inv(F);

    		//Correct the values on all meshes on a cell by cell basis

    		forAll(scalarVarsIncr[0], K)
    		{
    			scalarList newIncr(NT,0.0);

        		forAll(scalarVarsIncr,i)
        		{
        			forAll(scalarVarsIncr,j)
        			{
        				newIncr[i] += F[i][j]*scalarVarsIncr[j][K][celli];
        			}
        		}

        		forAll(scalarVarsIncr,i)
        		{
        			scalarVarsIncr[i][K][celli] = newIncr[i];
        		}
    		}

    		forAll(vectorVarsIncr[0], K)
    		{
    			List<vector> newIncr(NT, Zero);

        		forAll(vectorVarsIncr,i)
        		{
        			forAll(vectorVarsIncr,j)
        			{
        				newIncr[i] += F[i][j]*vectorVarsIncr[j][K][celli];
        			}
        		}

        		forAll(vectorVarsIncr,i)
        		{
        			vectorVarsIncr[i][K][celli] = newIncr[i];
        		}
    		}
    	}
    }
    else
    {
    	const labelList& cells = mesh_.cellZones()[cellZoneID_];

    	forAll (cells, index)
    	{
			label celli = cells[index];

    		// Compute the factorization matrix F

    		scalarList cellDeltaT(NT,0.0);
    		SquareMatrix<scalar> dTau(NT,0.0);

    		forAll(cellDeltaT, meshi)
    		{
    			cellDeltaT[meshi] = 1/rDeltaTField[meshi][celli];
				dTau[meshi][meshi] = cellDeltaT[meshi];
    		}

    		SquareMatrix<scalar> II(NT,0.0);
    		for (label l=0; l<NT; l++)
    		{
    			II[l][l] = 1.0;
    		}

    		SquareMatrix<scalar> f_0 = II + (dTau*this->D());
    		LUscalarMatrix f(f_0);
    		scalarSquareMatrix F(NT,0.0);
    		f.inv(F);

    		//Correct the values on all meshes on a cell by cell basis

    		forAll(scalarVarsIncr[0], K)
    		{
    			scalarList newIncr(NT,0.0);

        		forAll(scalarVarsIncr,i)
        		{
        			forAll(scalarVarsIncr,j)
        			{
        				newIncr[i] += F[i][j]*scalarVarsIncr[j][K][celli];
        			}
        		}

        		forAll(scalarVarsIncr,i)
        		{
        			scalarVarsIncr[i][K][celli] = newIncr[i];
        		}
    		}

    		forAll(vectorVarsIncr[0], K)
    		{
    			List<vector> newIncr(NT, Zero);

        		forAll(vectorVarsIncr,i)
        		{
        			forAll(vectorVarsIncr,j)
        			{
        				newIncr[i] += F[i][j]*vectorVarsIncr[j][K][celli];
        			}
        		}

        		forAll(vectorVarsIncr,i)
        		{
        			vectorVarsIncr[i][K][celli] = newIncr[i];
        		}
    		}
    	}
    }
}


void Foam::HBZone::reconstruct
(
	volVectorField& reconstrFld,
	PtrList<volVectorField>& fieldPtr,
	scalar actTime,
	bool tryCylCoords
) const
{
	const label nO = fieldPtr.size(); //number of snapshots
	const label nT = omegaList_.size(); //number of frequencies

	List<complex> E_1(nT);
	complex e(0,1);

	forAll (omegaList_, k)
	{
		e.Re() = Foam::cos(omegaList_[k]*actTime);
		e.Im() = Foam::sin(omegaList_[k]*actTime);
		E_1[k] = e;
	}

	complex t(0,0);
	List<complex> transfOp(nO, t);

	for (label m=0; m<nO; m++)
	{
		for (label n=0; n<nT; n++)
		{
			transfOp[m] += E_1[n]*E_[n][m];
		}
	}

	scalarList transfOpRe(nO, 0.0);

	forAll(transfOpRe, J)
	{
		transfOpRe[J] = transfOp[J].Re();
	}

    if (cellZoneID_ == -1)
    {
        return;
    }

    if (cylCoords_)
    {
    	dimensionedScalar smallRadius("smallRadius", dimLength, SMALL);

        if (cellZoneID_ == -2)
        {
        	label cellNo = mesh_.nCells();

        	vector axisHat = rotationAxis_/mag(rotationAxis_);

        	for (Foam::label celli=0; celli<cellNo; ++celli)
        	{
        		reconstrFld[celli] = Zero;
        		vector reconstrCyl(Zero);

        		for (label l=0; l<nO; l++)
        		{
					const fvMesh& meshL = fieldPtr[l].mesh();

					const vector& cellCentre = meshL.C()[celli];

					 // Radius vector in plane of rotation
					vector r(cellCentre - rotationCentre_);
					r -= (axisHat & r)*axisHat;
					const scalar magr(mag(r));
					const vector rHat(r/magr);

					scalar Ur = fieldPtr[l][celli] & rHat;
					scalar Uu = fieldPtr[l][celli] & (axisHat ^ rHat);
					scalar Uz = fieldPtr[l][celli] & axisHat;

        			vector UCyl(Ur, Uu, Uz);
        			reconstrCyl += transfOpRe[l]*UCyl;
        		}

				const vector& cellCentreAct = reconstrFld.mesh().C()[celli];

				// Radius vector in plane of rotation
				vector r(cellCentreAct - rotationCentre_);
				r -= (axisHat & r)*axisHat;
				const scalar magr(mag(r));
				const vector rHat(r/magr);

				vector reconstrCart = reconstrCyl.x()*rHat
								  + reconstrCyl.y()*(axisHat^rHat)
								  + reconstrCyl.z()*axisHat;

				reconstrFld[celli] = reconstrCart;
        	}
        }
        else
        {
        	const labelList& cells = mesh_.cellZones()[cellZoneID_];

        	vector axisHat = rotationAxis_/mag(rotationAxis_);

        	forAll (cells, i)
        	{
    			label celli = cells[i];

        		reconstrFld[celli] = Zero;
        		vector reconstrCyl(Zero);

				for (label l=0; l<nO; l++)
				{
					const fvMesh& meshL = fieldPtr[l].mesh();

					const vector& cellCentre = meshL.C()[celli];

					 // Radius vector in plane of rotation
					vector r(cellCentre - rotationCentre_);
					r -= (axisHat & r)*axisHat;
					const scalar magr(mag(r));
					const vector rHat(r/magr);

					scalar Ur = fieldPtr[l][celli] & rHat;
					scalar Uu = fieldPtr[l][celli] & (axisHat ^ rHat);
					scalar Uz = fieldPtr[l][celli] & axisHat;

					vector UCyl(Ur, Uu, Uz);
					reconstrCyl += transfOpRe[l]*UCyl;
				}

				const vector& cellCentreAct = reconstrFld.mesh().C()[celli];

				// Radius vector in plane of rotation
				vector r(cellCentreAct - rotationCentre_);
				r -= (axisHat & r)*axisHat;
				const scalar magr(mag(r));
				const vector rHat(r/magr);

				vector reconstrCart = reconstrCyl.x()*rHat
								  + reconstrCyl.y()*(axisHat^rHat)
								  + reconstrCyl.z()*axisHat;

				reconstrFld[celli] = reconstrCart;
        	}
        }
    }
    else
    {
    	reconstruct(reconstrFld, fieldPtr, actTime);
    }
}


void Foam::HBZone::writeData(Ostream& os) const
{
    os  << nl;
    os.beginBlock(name_);

    os.writeEntry("active", active_);
    os.writeEntry("cellZone", cellZoneName_);
    os.writeEntry("omegaList", omegaList_);
    os.writeEntry("cylCoords", cylCoords_);

    if (cylCoords_)
    {
    	os.writeEntry("rotationAxis", rotationAxis_);
    	os.writeEntry("rotationCentre", rotationCentre_);
    }

    os.writeEntry("EInv", EInv_);
    os.writeEntry("E", E_);
    os.writeEntry("D", D_);

    os.endBlock();
}


bool Foam::HBZone::read(const dictionary& dict)
{
    coeffs_ = dict;

    active_ = coeffs_.getOrDefault("active", true);
    coeffs_.readEntry("cellZone", cellZoneName_);
    cellZoneID_ = mesh_.cellZones().findZoneID(cellZoneName_);

    coeffs_.readEntry("cylCoords", cylCoords_);

    if (cylCoords_)
    {
    	coeffs_.readEntry("rotationAxis", rotationAxis_);
    	coeffs_.readEntry("rotationCentre", rotationCentre_);
    }

    return true;
}


// ************************************************************************* //
