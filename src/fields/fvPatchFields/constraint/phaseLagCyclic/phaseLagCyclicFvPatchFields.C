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

#include "phaseLagCyclicFvPatchFields.H"

#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

#include "psiThermo.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	makePatchFields(phaseLagCyclic);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::Field<scalar>> Foam::phaseLagCyclicFvPatchField<scalar>::phaseLaggedField() const
{
	const objectRegistry& allSubLevels = this->db().parent();
	const objectRegistry& subLevel0 = allSubLevels.lookupObject<objectRegistry>("subTimeLevel0");
	const HBZoneList& HB = subLevel0.lookupObject<IOHBZoneList>("HBProperties");

	scalar IBPA = 0.0;
	label HBZoneInstance = -1;

	if (this->phaseLagCyclicPatch().owner())
	{
		IBPA = this->phaseLagCyclicPatch().phaseLagCyclicPatch().IBPA();

		forAll (HB, i)
		{
			if (HB[i].name() == this->phaseLagCyclicPatch().phaseLagCyclicPatch().HBZoneName())
			{
				HBZoneInstance = i;
			}
		}
	}
	else
	{
		IBPA = -1*this->phaseLagCyclicPatch().phaseLagCyclicPatch().neighbPatch().IBPA();

		forAll (HB, i)
		{
			if (HB[i].name() == this->phaseLagCyclicPatch().phaseLagCyclicPatch().neighbPatch().HBZoneName())
			{
				HBZoneInstance = i;
			}
		}
	}

	tmp<Field<scalar>> tPhaseLag(new Field<scalar>(this->patchInternalField()));

	const dictionary& solDict = this->internalField().mesh().solutionDict();

	const label nT = solDict.subDict("harmonicBalance").getOrDefault<label>("instantsNumber",3);

	const RectangularMatrix<complex>& E = HB[HBZoneInstance].E();
	const RectangularMatrix<complex>& E_1 = HB[HBZoneInstance].EInv();

	const word& fieldName = this->internalField().name();

	const label& patchIndex = this->patch().index();
	const label& neighPatchIndex = this->phaseLagCyclicPatch().neighbPatchID();

	PtrList<Field<scalar>> perioFields(nT);

	forAll(perioFields,i)
	{
		word itrName = Foam::name(i);
		word timeLevel = "subTimeLevel" + itrName;

		const objectRegistry& subLeveli = allSubLevels.lookupObject<objectRegistry>(timeLevel);

		if (subLeveli.found(fieldName))
		{
			const volScalarField& subTimeField = subLeveli.lookupObject<volScalarField>(fieldName);

			if (!subTimeField.boundaryField().operator()(patchIndex)
				|| !subTimeField.boundaryField().operator()(neighPatchIndex))
			{
				Info<<"not correct"<<endl;
				return tPhaseLag;
			}

		    const scalarField& iField = subTimeField.primitiveField();
		    const labelUList& nbrFaceCells =
		        this->phaseLagCyclicPatch().neighbFvPatch().faceCells();

		    scalarField subTimeInternali(this->size());

	        forAll(*this, facei)
	        {
	            subTimeInternali[facei] = iField[nbrFaceCells[facei]];
	        }

			perioFields.set
			(
				i,
				new scalarField(subTimeInternali)
			);

		}
		else
		{
			if (!(subLeveli.found("U") && subLeveli.found("p") && subLeveli.found("rho")
				&& subLeveli.found("e") && subLeveli.found("thermophysicalProperties")))
			{
				return tPhaseLag;
			}
			else
			{
				const labelUList& nbrFaceCells =
						this->phaseLagCyclicPatch().neighbFvPatch().faceCells();

				const vectorField& Ui = subLeveli.lookupObject<volVectorField>("U").primitiveField();
				vectorField Up(Ui, nbrFaceCells);

				const scalarField& pi = subLeveli.lookupObject<volScalarField>("p").primitiveField();
				scalarField pp(pi, nbrFaceCells);

				const scalarField& rhoi = subLeveli.lookupObject<volScalarField>("rho").primitiveField();
				scalarField rhop(rhoi, nbrFaceCells);

				const scalarField& ei = subLeveli.lookupObject<volScalarField>("e").primitiveField();
				scalarField ep(ei, nbrFaceCells);

				const psiThermo& pThermo = subLeveli.lookupObject<psiThermo>("thermophysicalProperties");
				tmp< volScalarField > gammai = pThermo.gamma();
				scalarField	Gammai = gammai().primitiveField();
				scalarField gammap(Gammai, nbrFaceCells);

				tmp< volScalarField > psi = pThermo.psi();
				scalarField Psii = psi().primitiveField();
				scalarField psip(Psii, nbrFaceCells);


				if (fieldName == "(he+(0.5*magSqr(U)))")
				{
					scalarField subTimeInternali = ep + (0.5*magSqr(Up));

					perioFields.set
					(
						i,
						new scalarField(subTimeInternali)
					);
				}
				else if (fieldName == "((he+(0.5*magSqr(U)))+(p|rho))")
				{
					scalarField subTimeInternali = ep + (0.5*magSqr(Up)) + pp/rhop;

					perioFields.set
					(
						i,
						new scalarField(subTimeInternali)
					);
				}
				else if (fieldName == "U.component(0)")
				{
					scalarField subTimeInternali = Up.component(0);

					perioFields.set
					(
						i,
						new scalarField(subTimeInternali)
					);
				}
				else if (fieldName == "U.component(1)")
				{
					scalarField subTimeInternali = Up.component(1);

					perioFields.set
					(
						i,
						new scalarField(subTimeInternali)
					);

				}
				else if (fieldName == "U.component(2)")
				{
					scalarField subTimeInternali = Up.component(2);

					perioFields.set
					(
						i,
						new scalarField(subTimeInternali)
					);
				}
				else if (fieldName == "gamma")
				{
					scalarField subTimeInternali(gammap);

					perioFields.set
					(
						i,
						new scalarField(subTimeInternali)
					);
				}
				else if (fieldName == "c")
				{
					scalarField subTimeInternali = sqrt(gammap/psip);

					perioFields.set
					(
						i,
						new scalarField(subTimeInternali)
					);
				}
				else if (fieldName == "sqrt((gamma|thermo:psi))")
				{
					scalarField subTimeInternali = sqrt(gammap/psip);

					perioFields.set
					(
						i,
						new scalarField(subTimeInternali)
					);
				}
				else if (fieldName == "H")
				{
					scalarField subTimeInternali = ep + (0.5*magSqr(Up)) + pp/rhop;

					perioFields.set
					(
						i,
						new scalarField(subTimeInternali)
					);
				}
				else
				{
					Info<<"not found field Name "<<fieldName<<" for time level "<<timeLevel<<endl;

					return tPhaseLag;
				}
			}
		}
	}

	label nF = E.m();
	label nH = (nF-1)/2;

	complex t(0,0);
	SquareMatrix<complex> M(nF, t);
	SquareMatrix<complex> d(nT, t);
	RectangularMatrix<complex> temp0(nF, nT, t);

	M[0][0] = complex(1,0);

	for (int n = 1; n <= nH; n++)
	{
		t.Re() = Foam::cos(n*IBPA);
		t.Im() = Foam::sin(n*IBPA);
		M[n][n] = t;
		M[nF-n][nF-n] = t.conjugate();
	}

	for (label l=0; l<nF; l++)
	{
		for (label m=0; m<nT; m++)
		{
			for (label n=0; n<nF; n++)
			{
				temp0[l][m] += M[l][n]*E[n][m];
			}
		}
	}

	for (label l=0; l<nT; l++)
	{
		for (label m=0; m<nT; m++)
		{
			for (label n=0; n<nF; n++)
			{
				d[l][m] += E_1[l][n]*temp0[n][m];
			}
		}
	}

	SquareMatrix<scalar> D(nT, 0.0);
	const Identity<scalar> ii;
	SquareMatrix<scalar> Id(nT, ii);

	for (int i = 0; i < nT; i++)
	{
		for (int j = 0; j < nT; j++)
		{
			D[i][j] = d[i][j].Re();
		}
	}

	Field<scalar>& phaseLag = tPhaseLag.ref();

	bool cylCoords = false;

	if (this->phaseLagCyclicPatch().owner())
	{
		cylCoords = this->phaseLagCyclicPatch().phaseLagCyclicPatch().cylCoords();
	}
	else
	{
		cylCoords = this->phaseLagCyclicPatch().phaseLagCyclicPatch().neighbPatch().cylCoords();
	}

	if (cylCoords)
	{
		if (fieldName.startsWith("U.component"))
		{
			dimensionedScalar smallRadius("smallRadius", dimLength, SMALL);
			phaseLag = Zero;

			PtrList<vectorField> UPtr(nT);

			forAll(perioFields, J)
			{
				word itrName = Foam::name(J);
				word timeLevel = "subTimeLevel" + itrName;

				const objectRegistry& subLeveli = allSubLevels.lookupObject<objectRegistry>(timeLevel);

				const labelUList& nbrFaceCells = this->phaseLagCyclicPatch().neighbFvPatch().faceCells();
				const volVectorField& Ui = subLeveli.lookupObject<volVectorField>("U");

				vectorField UPatch(Ui.primitiveField(), nbrFaceCells);

				UPtr.set(J, new vectorField(UPatch));
			}

			const vector& rotAxis = this->phaseLagCyclicPatch().phaseLagCyclicPatch().rotationAxis();
		    const vector axisHat = rotAxis/mag(rotAxis);

		    const point& origin = this->phaseLagCyclicPatch().phaseLagCyclicPatch().rotationCentre();


			forAll(phaseLag, facei)
			{
				vector sourceCyl(Zero);

				forAll(perioFields, J)
				{
					word itrName = Foam::name(J);
					word timeLevel = "subTimeLevel" + itrName;

					const objectRegistry& subLeveli = allSubLevels.lookupObject<objectRegistry>(timeLevel);

					const volVectorField& Ui = subLeveli.lookupObject<volVectorField>("U");

					const vectorField& faceCenters = Ui.mesh().boundary()[neighPatchIndex].Cf();

				    // Radius vector in plane of rotation
				    vector r(faceCenters[facei] - origin);
				    r -= (axisHat & r)*axisHat;
				    const scalar magr(mag(r));
				    const vector rHat(r/magr);

					const vectorField& Up = UPtr[J];

					scalar Ur = Up[facei] & rHat;
					scalar Uu = Up[facei] & (axisHat ^ rHat);
					scalar Uz = Up[facei] & axisHat;

					vector UCyl(Ur, Uu, Uz);

					sourceCyl += D[subTimeLevel_][J]*UCyl;
				}

				const vectorField& faceCentersAct = phaseLagCyclicPatch().neighbPatch().Cf();

				// Radius vector in plane of rotation
				vector r(faceCentersAct[facei] - origin);
				r -= (axisHat & r)*axisHat;
				const scalar magr(mag(r));
				const vector rHat(r/magr);

			    vector sourceCart = sourceCyl.x()*rHat
			    				  + sourceCyl.y()*(axisHat^rHat)
								  + sourceCyl.z()*axisHat;

			    if (doTransform())
			    {
			    	sourceCart = transform(forwardT()[0], sourceCart);
			    }

			    if (fieldName == "U.component(0)")
			    {
			    	 phaseLag[facei] += sourceCart.x();
			    }
			    else if (fieldName == "U.component(1)")
			    {
			    	phaseLag[facei] += sourceCart.y();
			    }
			    else
				{
					phaseLag[facei] += sourceCart.z();
				}
			}
		}
		else
		{
			phaseLag = Zero;

			forAll(phaseLag, facei)
			{
				forAll(perioFields, J)
				{
					phaseLag[facei] += D[subTimeLevel_][J]*perioFields[J][facei];
				}
			}
		}
	}
	else
	{
		phaseLag = Zero;

		forAll(phaseLag, facei)
		{
			forAll(perioFields, J)
			{
				phaseLag[facei] += D[subTimeLevel_][J]*perioFields[J][facei];
			}
		}
	}

	return tPhaseLag;
}


template<>
Foam::tmp<Foam::vectorField> Foam::phaseLagCyclicFvPatchField<vector>::phaseLaggedField() const
{
	const objectRegistry& allSubLevels = this->db().parent();
  	const objectRegistry& subLevel0 = allSubLevels.lookupObject<objectRegistry>("subTimeLevel0");
	const HBZoneList& HB = subLevel0.lookupObject<IOHBZoneList>("HBProperties");

	scalar IBPA = 0.0;
	label HBZoneInstance = -1;

	if (this->phaseLagCyclicPatch().owner())
	{
		IBPA = this->phaseLagCyclicPatch().phaseLagCyclicPatch().IBPA();

		forAll (HB, i)
		{
			if (HB[i].name() == this->phaseLagCyclicPatch().phaseLagCyclicPatch().HBZoneName())
			{
				HBZoneInstance = i;
			}
		}
	}
	else
	{
		IBPA = -1*this->phaseLagCyclicPatch().phaseLagCyclicPatch().neighbPatch().IBPA();

		forAll (HB, i)
		{
			if (HB[i].name() == this->phaseLagCyclicPatch().phaseLagCyclicPatch().neighbPatch().HBZoneName())
			{
				HBZoneInstance = i;
			}
		}
	}

	if (HBZoneInstance == -1)
	{
        FatalErrorInFunction
            << " Specified HBZoneName not found "
            << exit(FatalIOError);
	}

	tmp<vectorField> tPhaseLag(new vectorField(this->patchInternalField()));

	const dictionary& solDict = this->internalField().mesh().solutionDict();

	const label nT = solDict.subDict("harmonicBalance").getOrDefault<label>("instantsNumber",3);

	const RectangularMatrix<complex>& E = HB[HBZoneInstance].E();
	const RectangularMatrix<complex>& E_1 = HB[HBZoneInstance].EInv();

	const word& fieldName = this->internalField().name();

	const label& patchIndex = this->patch().index();
	const label& neighPatchIndex = this->phaseLagCyclicPatch().neighbPatchID();

	PtrList<vectorField> perioFields(nT);

	forAll(perioFields,i)
	{
		word itrName = Foam::name(i);
		word timeLevel = "subTimeLevel" + itrName;

		const objectRegistry& subLeveli = allSubLevels.lookupObject<objectRegistry>(timeLevel);

		if (subLeveli.found(fieldName))
		{
			const volVectorField& subTimeField =
					subLeveli.lookupObject<volVectorField>(fieldName);

			if (!subTimeField.boundaryField().operator()(patchIndex)
				|| !subTimeField.boundaryField().operator()(neighPatchIndex))
			{
				Info<<"not correct"<<endl;
				return tPhaseLag;
			}

		    const vectorField& iField = subTimeField.primitiveField();
		    const labelUList& nbrFaceCells =
		        this->phaseLagCyclicPatch().neighbFvPatch().faceCells();

		    vectorField subTimeInternali(this->size());

	        forAll(*this, facei)
	        {
	            subTimeInternali[facei] = iField[nbrFaceCells[facei]];
	        }

			perioFields.set
			(
				i,
				new vectorField(subTimeInternali)
			);
		}
		else
		{
			Info<<"not correct"<<endl;
			return tPhaseLag;
		}
	}

	label nF = E.m();
	label nH = (nF-1)/2;

	complex t(0,0);
	SquareMatrix<complex> M(nF, t);
	SquareMatrix<complex> d(nT, t);
	RectangularMatrix<complex> temp0(nF, nT, t);

	M[0][0] = complex(1,0);

	for (int n = 1; n <= nH; n++)
	{
		t.Re() = Foam::cos(n*IBPA);
		t.Im() = Foam::sin(n*IBPA);
		M[n][n] = t;
		M[nF-n][nF-n] = t.conjugate();
	}

	for (label l=0; l<nF; l++)
	{
		for (label m=0; m<nT; m++)
		{
			for (label n=0; n<nF; n++)
			{
				temp0[l][m] += M[l][n]*E[n][m];
			}
		}
	}

	for (label l=0; l<nT; l++)
	{
		for (label m=0; m<nT; m++)
		{
			for (label n=0; n<nF; n++)
			{
				d[l][m] += E_1[l][n]*temp0[n][m];
			}
		}
	}

	SquareMatrix<scalar> D(nT, 0.0);
	const Identity<scalar> ii;
	SquareMatrix<scalar> Id(nT, ii);

	for (int i = 0; i < nT; i++)
	{
		for (int j = 0; j < nT; j++)
		{
			D[i][j] = d[i][j].Re();
		}
	}

	vectorField& phaseLag = tPhaseLag.ref();

	phaseLag = Zero;

	bool cylCoords = false;

	if (this->phaseLagCyclicPatch().owner())
	{
		cylCoords = this->phaseLagCyclicPatch().phaseLagCyclicPatch().cylCoords();
	}
	else
	{
		cylCoords = this->phaseLagCyclicPatch().phaseLagCyclicPatch().neighbPatch().cylCoords();
	}

	if (cylCoords)
	{
		dimensionedScalar smallRadius("smallRadius", dimLength, SMALL);

		const vector& rotAxis = this->phaseLagCyclicPatch().phaseLagCyclicPatch().rotationAxis();
		const vector axisHat = rotAxis/mag(rotAxis);
		const point& origin = this->phaseLagCyclicPatch().phaseLagCyclicPatch().rotationCentre();

		forAll(phaseLag, facei)
		{
			vector sourceCyl(Zero);

			forAll(perioFields, J)
			{
				word itrName = Foam::name(J);
				word timeLevel = "subTimeLevel" + itrName;

				const objectRegistry& subLeveli = allSubLevels.lookupObject<objectRegistry>(timeLevel);

				const volVectorField& Ui = subLeveli.lookupObject<volVectorField>("U");

				const vectorField& faceCenters = Ui.mesh().boundary()[neighPatchIndex].Cf();

				// Radius vector in plane of rotation
				vector r(faceCenters[facei] - origin);
				r -= (axisHat & r)*axisHat;
				const scalar magr(mag(r));
				const vector rHat(r/magr);

				const vectorField& Up = perioFields[J];

				scalar Ur = Up[facei] & rHat;
				scalar Uu = Up[facei] & (axisHat ^ rHat);
				scalar Uz = Up[facei] & axisHat;

				vector UCyl(Ur, Uu, Uz);

				sourceCyl += D[subTimeLevel_][J]*UCyl;
			}

			const vectorField& faceCentersAct = phaseLagCyclicPatch().neighbPatch().Cf();

			// Radius vector in plane of rotation
			vector r(faceCentersAct[facei] - origin);
			r -= (axisHat & r)*axisHat;
			const scalar magr(mag(r));
			const vector rHat(r/magr);

			vector sourceCart = sourceCyl.x()*rHat
							  + sourceCyl.y()*(axisHat^rHat)
							  + sourceCyl.z()*axisHat;

			if (doTransform())
			{
				sourceCart = transform(forwardT()[0], sourceCart);
			}

			phaseLag[facei] = sourceCart;
		}
	}
	else
	{
		forAll(phaseLag, facei)
		{
			forAll(perioFields, J)
			{
				phaseLag[facei] += D[subTimeLevel_][J]*perioFields[J][facei];
			}
		}

	    if (this->doTransform())
	    {
	    	forAll(phaseLag, facei)
			{
	    		phaseLag[facei] = transform
	    		(
	    			this->forwardT()[0], phaseLag[facei]
	    		);
			}
	    }
	}

	return tPhaseLag;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
