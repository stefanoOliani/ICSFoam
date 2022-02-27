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

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::HBZone::addSource
(
	PtrList<Field<Type>>& source,
	PtrList<PtrList<GeometricField<Type, fvPatchField, volMesh>>>& vars,
	label varNum
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const scalarField& V = mesh_.V();

    if (cellZoneID_ == -2)
    {
    	label cellNo = mesh_.nCells();

    	forAll(source, J)
    	{
    		source[J] = Zero;
    	}

    	for (Foam::label celli=0; celli<cellNo; ++celli)
    	{
    		forAll (source, J)
    		{
    			forAll (source, K)
				{
					source[J][celli] -= V[celli]*this->D()[J][K]*vars[K][varNum][celli];
				}
    		}
    	}
    }
    else
    {
    	const labelList& cells = mesh_.cellZones()[cellZoneID_];

    	forAll (cells, i)
    	{
			label celli = cells[i];

    		forAll (source, J)
    		{
    			source[J][celli] = Zero;

    			forAll (source, K)
				{
					source[J][celli] -= V[celli]*this->D()[J][K]*vars[K][varNum][celli];
				}
    		}
    	}
    }
}

template<class Type>
void Foam::HBZone::addSource
(
	fvMatrix<Type>& eqn,
	PtrList<GeometricField<Type, fvPatchField, volMesh>>& fieldPtr,
	label instantNo
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const scalarField& V = mesh_.V();
    Field<Type>& Source = eqn.source();

    if (cellZoneID_ == -2)
    {
    	label cellNo = mesh_.nCells();

    	for (Foam::label celli=0; celli<cellNo; ++celli)
    	{
    		forAll (fieldPtr, J)
    		{
    			Source[celli] -= V[celli]*this->D()[instantNo][J]*fieldPtr[J][celli];
    		}
    	}
    }
    else
    {
    	const labelList& cells = mesh_.cellZones()[cellZoneID_];

    	forAll (cells, i)
    	{
			label celli = cells[i];

    		forAll (fieldPtr, J)
    		{
    			Source[celli] -= V[celli]*this->D()[instantNo][J]*fieldPtr[J][celli];
    		}
    	}
    }
}

template<class Type>
void Foam::HBZone::addSource
(
	const volScalarField& rho,
	fvMatrix<Type>& eqn,
	PtrList<GeometricField<Type, fvPatchField, volMesh>>& fieldPtr,
	label instantNo
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const scalarField& V = mesh_.V();
    Field<Type>& Source = eqn.source();

    if (cellZoneID_ == -2)
    {
    	label cellNo = mesh_.nCells();

    	for (Foam::label celli=0; celli<cellNo; ++celli)
    	{
    		forAll (fieldPtr, J)
    		{
    			Source[celli] -= V[celli]*rho[celli]*this->D()[instantNo][J]*fieldPtr[J];
    		}
    	}
    }
    else
    {
    	const labelList& cells = mesh_.cellZones()[cellZoneID_];

    	forAll (cells, i)
    	{
			label celli = cells[i];

    		forAll (fieldPtr, J)
    		{
    			Source[celli] -= V[celli]*rho[celli]*this->D()[instantNo][J]*fieldPtr[J];
    		}
    	}
    }
}


template<class Type>
void Foam::HBZone::factorizationStep
(
	PtrList<GeometricField<Type, fvPatchField, volMesh>>& fieldPtr,
	PtrList<scalarField>& deltaTField
) const
{
	label NT = fieldPtr.size();

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
    			cellDeltaT[meshi] = deltaTField[meshi][celli];
				dTau[meshi][meshi] = cellDeltaT[meshi];
    		}

    		SquareMatrix<scalar> II(NT,0.0);
    		for (label l=0; l<NT; l++)
    		{
    			II[l][l] = 1.0;
    		}

    		SquareMatrix<scalar> f_0 = II + dTau*this->D();
    		LUscalarMatrix f(f_0);
    		scalarSquareMatrix F(NT,0.0);
    		f.inv(F);

    		//Correct the values on all meshes on a cell by cell basis

    		List<Type> T_cellNew(NT,Zero);
    		List<Type> T_cell(NT,Zero);
    		List<Type> T_cellOld(NT,Zero);

    		forAll(T_cell, meshi)
    		{
    			T_cell[meshi] = fieldPtr[meshi][celli];
    			T_cellOld[meshi] = fieldPtr[meshi].oldTime()[celli];
    		}

    		List<Type> T_toCorrect(NT, Zero);

    		T_toCorrect = T_cell - T_cellOld;

    		forAll(T_cellNew,i)
    		{
    			forAll(T_cellNew,j)
    			{
    				T_cellNew[i] += F[i][j]*T_toCorrect[j];
    			}
    		}

    		T_cellNew = T_cellNew + T_cellOld;

    		// Update values on each mesh
    		forAll(fieldPtr, meshi)
    		{
    			fieldPtr[meshi][celli] = T_cellNew[meshi];
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
    			cellDeltaT[meshi] = deltaTField[meshi][celli];
				dTau[meshi][meshi] = cellDeltaT[meshi];
    		}

    		SquareMatrix<scalar> II(NT,0.0);
    		for (label l=0; l<NT; l++)
    		{
    			II[l][l] = 1.0;
    		}

    		SquareMatrix<scalar> f_0 = II + dTau*this->D();
    		LUscalarMatrix f(f_0);
    		scalarSquareMatrix F(NT,0.0);
    		f.inv(F);

    		//Correct the values on all meshes on a cell by cell basis

    		List<Type> T_cellNew(NT,Zero);
    		List<Type> T_cell(NT,Zero);
    		List<Type> T_cellOld(NT,Zero);

    		forAll(T_cell, meshi)
    		{
    			T_cell[meshi] = fieldPtr[meshi][celli];
    			T_cellOld[meshi] = fieldPtr[meshi].oldTime()[celli];
    		}

    		List<Type> T_toCorrect(NT, Zero);

    		T_toCorrect = T_cell - T_cellOld;

    		forAll(T_cellNew,i)
    		{
    			forAll(T_cellNew,j)
    			{
    				T_cellNew[i] += F[i][j]*T_toCorrect[j];
    			}
    		}

    		T_cellNew = T_cellNew + T_cellOld;

    		// Update values on each mesh
    		forAll(fieldPtr, meshi)
    		{
    			fieldPtr[meshi][celli] = T_cellNew[meshi];
    		}

    	}
    }
}


template<class Type>
void Foam::HBZone::reconstruct
(
	GeometricField<Type, fvPatchField, volMesh>& reconstrFld,
	PtrList<GeometricField<Type, fvPatchField, volMesh>>& fieldPtr,
	scalar actTime
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


    if (cellZoneID_ == -1)
    {
        return;
    }

    if (cellZoneID_ == -2)
    {
    	const labelList cells(mesh_.nCells());

    	label cellNo = mesh_.nCells();

    	for (Foam::label celli=0; celli<cellNo; ++celli)
    	{
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

    		reconstrFld[celli] = Zero;

    		for (label l=0; l<nO; l++)
    		{
    			reconstrFld[celli] += transfOpRe[l]*fieldPtr[l][celli];
    		}
    	}
    }
    else
    {
    	const labelList& cells = mesh_.cellZones()[cellZoneID_];

    	forAll (cells, i)
    	{
			label celli = cells[i];

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

    		reconstrFld[celli] = Zero;

    		for (label l=0; l<nO; l++)
    		{
    			reconstrFld[celli] += transfOpRe[l]*fieldPtr[l][celli];
    		}
    	}
    }

}


// ************************************************************************* //
