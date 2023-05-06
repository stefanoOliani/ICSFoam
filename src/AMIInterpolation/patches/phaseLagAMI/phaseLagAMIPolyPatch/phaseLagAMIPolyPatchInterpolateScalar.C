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

#include "phaseLagAMIPolyPatch.H"
#include "psiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::scalarField >
Foam::phaseLagAMIPolyPatch::expandData(const scalarField& pf, const word& fieldName) const
{
    // Check and expand the field from patch size to zone size
    if (pf.size() != this->size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > phaseLagAMIPolyPatch::expandData"
            "("
            "    const Field<Type>& pf"
            ") const"
        )   << "Incorrect patch field size.  Field size: "
            << pf.size() << " patch size: " << this->size()
            << abort(FatalError);
    }

	const objectRegistry& allSubLevels = this->boundaryMesh().mesh().objectRegistry::parent();
  	const objectRegistry& subLevel0 = allSubLevels.lookupObject<objectRegistry>("subTimeLevel0");
	const HBZoneList& HB = subLevel0.lookupObject<IOHBZoneList>("HBProperties");

	label HBZoneInstance = -1;

	forAll (HB, i)
	{
		if (HB[i].name() == HBZoneName_)
		{
			HBZoneInstance = i;
		}
	}

	const label nT = HB.selectedSnapshots().size();

	const RectangularMatrix<complex>& E = HB[HBZoneInstance].E();
	const RectangularMatrix<complex>& E_1 = HB[HBZoneInstance].EInv();

	const label& patchIndex = this->index();
	const label& neighPatchIndex = this->neighbPatchID();

	PtrList<scalarField> perioFields(nT);

    tmp<scalarField> texpandField
    (
        new scalarField(pf)
    );

	forAll(perioFields,i)
	{
		word itrName = Foam::name(i);
		word timeLevel = "subTimeLevel" + itrName;

		const objectRegistry& subLeveli = allSubLevels.lookupObject<objectRegistry>(timeLevel);

		if(subLeveli.found(fieldName))
		{
			const volScalarField& subTimeField =
					subLeveli.lookupObject<volScalarField>(fieldName);

			if (!subTimeField.boundaryField().operator()(patchIndex)
				|| !subTimeField.boundaryField().operator()(neighPatchIndex))
			{
				const scalarField subTimeInternali = pf;

				perioFields.set
				(
					i,
					new scalarField(subTimeInternali)
				);
			}
			else
			{
				const scalarField subTimeInternali
							= subTimeField.boundaryField()[patchIndex].patchInternalField();

				perioFields.set
				(
					i,
					new scalarField(subTimeInternali)
				);
			}
		}
		else
		{
			if (!(subLeveli.found("U") && subLeveli.found("p") && subLeveli.found("rho")
					&& subLeveli.found("e") && subLeveli.found("thermophysicalProperties")))
			{
				const scalarField subTimeInternali = pf;

				perioFields.set
				(
					i,
					new scalarField(subTimeInternali)
				);
			}
			else
			{
				const labelUList& faceCells = this->faceCells();

				const vectorField& Ui = subLeveli.lookupObject<volVectorField>("U").primitiveField();
				vectorField Up(Ui, faceCells);

				const scalarField& pi = subLeveli.lookupObject<volScalarField>("p").primitiveField();
				scalarField pp(pi, faceCells);

				const scalarField& rhoi = subLeveli.lookupObject<volScalarField>("rho").primitiveField();
				scalarField rhop(rhoi, faceCells);

				const scalarField& ei = subLeveli.lookupObject<volScalarField>("e").primitiveField();
				scalarField ep(ei, faceCells);

				const psiThermo& pThermo = subLeveli.lookupObject<psiThermo>("thermophysicalProperties");
				tmp< volScalarField > gammai = pThermo.gamma();
				scalarField	Gammai = gammai().primitiveField();
				scalarField gammap(Gammai, faceCells);

				tmp< volScalarField > psi = pThermo.psi();
				scalarField Psii = psi().primitiveField();
				scalarField psip(Psii, faceCells);

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
					scalarField subTimeInternali = gammap;

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
					scalarField subTimeInternali = pf;

					perioFields.set
					(
						i,
						new scalarField(subTimeInternali)
					);
				}
			}
		}
	}

    // Get the periodic patch
    const coupledPolyPatch& periodicPatch
    (
        refCast<const coupledPolyPatch>
        (
            boundaryMesh()[periodicPatchID()]
        )
    );

    scalarField& expandField = texpandField.ref();
	scalarField transfField(pf.size(), Zero);

	label nF = E.m();
	label nH = (nF-1)/2;

    for (label copyI = 0; copyI < nTransformsBwd_; copyI++)
    {
		complex t(0,0);
		SquareMatrix<complex> M(nF, t);
		SquareMatrix<complex> d(nT, t);
		RectangularMatrix<complex> temp0(nF, nT, t);

		M[0][0] = complex(1,0);

		for (int n = 1; n <= nH; n++)
		{
			t.Re() = Foam::cos(-(copyI+1)*n*IBPA_);
			t.Im() = Foam::sin(-(copyI+1)*n*IBPA_);
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

		if (cylCoords_)
		{
			if (fieldName.startsWith("U.component"))
			{
				dimensionedScalar smallRadius("smallRadius", dimLength, SMALL);

				PtrList<vectorField> UPtr(nT);

				forAll(perioFields, J)
				{
					word itrName = Foam::name(J);
					word timeLevel = "subTimeLevel" + itrName;

					const objectRegistry& subLeveli = allSubLevels.lookupObject<objectRegistry>(timeLevel);

					const labelUList& faceCells = this->faceCells();
					const volVectorField& Ui = subLeveli.lookupObject<volVectorField>("U");
					vectorField UPatch(Ui.primitiveField(), faceCells);

					UPtr.set(J, new vectorField(UPatch));
				}

				const vector& rotAxis = HB[HBZoneInstance].rotationAxis();
			    const vector axisHat = rotAxis/mag(rotAxis);

			    const point& origin = HB[HBZoneInstance].rotationCentre();

				forAll(pf, facei)
				{
					vector sourceCyl(Zero);

					forAll(perioFields, J)
					{
						word itrName = Foam::name(J);
						word timeLevel = "subTimeLevel" + itrName;

						const objectRegistry& subLeveli = allSubLevels.lookupObject<objectRegistry>(timeLevel);

						const volVectorField& Ui = subLeveli.lookupObject<volVectorField>("U");

						const vectorField& Up = UPtr[J];

						const vectorField& faceCenters = Ui.mesh().boundary()[patchIndex].Cf();

						// Radius vector in plane of rotation
						vector r(faceCenters[facei] - origin);
						r -= (axisHat & r)*axisHat;
						const scalar magr(mag(r));
						const vector rHat(r/magr);

						scalar Ur = Up[facei] & rHat;
						scalar Uu = Up[facei] & (axisHat ^ rHat);
						scalar Uz = Up[facei] & axisHat;

						vector UCyl(Ur, Uu, Uz);

						sourceCyl += D[subTimeLevel_][J]*UCyl;
					}

					const vectorField& faceCentersAct = this->faceCentres();

					// Radius vector in plane of rotation
					vector r(faceCentersAct[facei] - origin);
					r -= (axisHat & r)*axisHat;
					const scalar magr(mag(r));
					const vector rHat(r/magr);

					vector sourceCart = sourceCyl.x()*rHat
									  + sourceCyl.y()*(axisHat^rHat)
									  + sourceCyl.z()*axisHat;

			        // Calculate transform
					tensor curTransform = periodicPatch.reverseT()[0];

					for (label instI = 0; instI < copyI; instI++)
					{
						curTransform &= curTransform;
					}

			    	if (curTransform.size())
			    	{
			    		 sourceCart = Foam::transform(curTransform, sourceCart);
			    	}

				    if (fieldName == "U.component(0)")
				    {
				    	transfField[facei] += sourceCart.x();
				    }
				    else if (fieldName == "U.component(1)")
				    {
				    	transfField[facei] += sourceCart.y();
				    }
				    else
				    {
				    	transfField[facei] += sourceCart.z();
				    }
				}
			}
			else
			{
				forAll(pf, facei)
				{
					forAll(perioFields, J)
					{
						transfField[facei] += D[subTimeLevel_][J]*perioFields[J][facei];
					}
				}
			}
		}
		else
		{
			forAll(pf, facei)
			{
				forAll(perioFields, J)
				{
					transfField[facei] += D[subTimeLevel_][J]*perioFields[J][facei];
				}
			}
		}

		expandField.append(transfField);

	    transfField = Zero;
    }

    for (label copyI = 0; copyI < nTransformsFwd_; copyI++)
    {
		complex t(0,0);
		SquareMatrix<complex> M(nF, t);
		SquareMatrix<complex> d(nT, t);
		RectangularMatrix<complex> temp0(nF, nT, t);

		M[0][0] = complex(1,0);

		for (int n = 1; n <= nH; n++)
		{
			t.Re() = Foam::cos((copyI+1)*n*IBPA_);
			t.Im() = Foam::sin((copyI+1)*n*IBPA_);
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

		if (cylCoords_)
		{
			if (fieldName.startsWith("U.component"))
			{
				dimensionedScalar smallRadius("smallRadius", dimLength, SMALL);

				PtrList<vectorField> UPtr(nT);

				forAll(perioFields, J)
				{
					word itrName = Foam::name(J);
					word timeLevel = "subTimeLevel" + itrName;

					const objectRegistry& subLeveli = allSubLevels.lookupObject<objectRegistry>(timeLevel);

					const labelUList& faceCells = this->faceCells();
					const volVectorField& Ui = subLeveli.lookupObject<volVectorField>("U");
					vectorField UPatch(Ui.primitiveField(), faceCells);

					UPtr.set(J, new vectorField(UPatch));
				}

				const vector& rotAxis = HB[HBZoneInstance].rotationAxis();
			    const vector axisHat = rotAxis/mag(rotAxis);

			    const point& origin = HB[HBZoneInstance].rotationCentre();

				forAll(pf, facei)
				{
					vector sourceCyl(Zero);

					forAll(perioFields, J)
					{
						word itrName = Foam::name(J);
						word timeLevel = "subTimeLevel" + itrName;

						const objectRegistry& subLeveli = allSubLevels.lookupObject<objectRegistry>(timeLevel);

						const volVectorField& Ui = subLeveli.lookupObject<volVectorField>("U");

						const vectorField& Up = UPtr[J];

						const vectorField& faceCenters = Ui.mesh().boundary()[patchIndex].Cf();

						// Radius vector in plane of rotation
						vector r(faceCenters[facei] - origin);
						r -= (axisHat & r)*axisHat;
						const scalar magr(mag(r));
						const vector rHat(r/magr);

						scalar Ur = Up[facei] & rHat;
						scalar Uu = Up[facei] & (axisHat ^ rHat);
						scalar Uz = Up[facei] & axisHat;

						vector UCyl(Ur, Uu, Uz);

						sourceCyl += D[subTimeLevel_][J]*UCyl;
					}

					const vectorField& faceCentersAct = this->faceCentres();

					// Radius vector in plane of rotation
					vector r(faceCentersAct[facei] - origin);
					r -= (axisHat & r)*axisHat;
					const scalar magr(mag(r));
					const vector rHat(r/magr);

					vector sourceCart = sourceCyl.x()*rHat
									  + sourceCyl.y()*(axisHat^rHat)
									  + sourceCyl.z()*axisHat;

			        // Calculate transform
					tensor curTransform = periodicPatch.forwardT()[0];

					for (label instI = 0; instI < copyI; instI++)
					{
						curTransform &= curTransform;
					}

			    	if (curTransform.size())
			    	{
			    		 sourceCart = Foam::transform(curTransform, sourceCart);
			    	}

				    if (fieldName == "U.component(0)")
				    {
				    	transfField[facei] += sourceCart.x();
				    }
				    else if (fieldName == "U.component(1)")
				    {
				    	transfField[facei] += sourceCart.y();
				    }
				    else
				    {
				    	transfField[facei] += sourceCart.z();
				    }
				}
			}
			else
			{
				forAll(pf, facei)
				{
					forAll(perioFields, J)
					{
						transfField[facei] += D[subTimeLevel_][J]*perioFields[J][facei];
					}
				}
			}
		}
		else
		{
			forAll(pf, facei)
			{
				forAll(perioFields, J)
				{
					transfField[facei] += D[subTimeLevel_][J]*perioFields[J][facei];
				}
			}
		}

		 expandField.append(transfField);

		 transfField = Zero;
    }

    return texpandField;
}


// ************************************************************************* //
