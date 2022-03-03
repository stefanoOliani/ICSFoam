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

#include "phaseLagCyclicFvPatchField.H"
#include "IOHBZoneList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::phaseLagCyclicFvPatchField<Type>::phaseLagCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
	cyclicFvPatchField<Type>(p, iF),
	subTimeLevel_(0),
	phaseLagCyclicPatch_(refCast<const phaseLagCyclicFvPatch>(p))
{
	const objectRegistry& allSubLevels = this->db().parent();

	if (allSubLevels.found("subTimeLevel0"))
	{
	  	const fvMesh& subLevel0 = allSubLevels.lookupObject<fvMesh>("subTimeLevel0");
	    const dictionary& solDict = subLevel0.solutionDict();

	    if (solDict.findDict("harmonicBalance"))
	    {
	        const scalar NT = solDict.subDict("harmonicBalance").getOrDefault<scalar>("instantsNumber",3);

	        for (label i=0; i<NT; i++)
	        {
	        	word itrName = Foam::name(i);
	        	word timeLevel = "subTimeLevel" + itrName;

	        	if (timeLevel == this->db().name())
	        	{
	        		subTimeLevel_ = i;
	        		break;
	        	}
	        }
	    }
	}
}


template<class Type>
Foam::phaseLagCyclicFvPatchField<Type>::phaseLagCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
	cyclicFvPatchField<Type>(p, iF),
	subTimeLevel_(0),
	phaseLagCyclicPatch_(refCast<const phaseLagCyclicFvPatch>(p, dict))

{
    if (!isA<phaseLagCyclicFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        this->evaluate(Pstream::commsTypes::blocking);
    }

    const objectRegistry& allSubLevels = this->db().parent();

    if (allSubLevels.foundObject<dictionary>("fvSolution"))
    {
        const dictionary& solDict = allSubLevels.lookupObject<dictionary>("fvSolution");

        const scalar NT = solDict.subDict("harmonicBalance").getOrDefault<scalar>("instantsNumber",3);

        for (label i=0; i<NT; i++)
        {
        	word itrName = Foam::name(i);
        	word timeLevel = "subTimeLevel" + itrName;

        	if (timeLevel == this->db().name())
        	{
        		subTimeLevel_ = i;
        		break;
        	}
        }
    }
}


template<class Type>
Foam::phaseLagCyclicFvPatchField<Type>::phaseLagCyclicFvPatchField
(
    const phaseLagCyclicFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
	cyclicFvPatchField<Type>(ptf, p, iF, mapper),
	subTimeLevel_(ptf.subTimeLevel_),
	phaseLagCyclicPatch_(refCast<const phaseLagCyclicFvPatch>(p))
{}


template<class Type>
Foam::phaseLagCyclicFvPatchField<Type>::phaseLagCyclicFvPatchField
(
    const phaseLagCyclicFvPatchField<Type>& ptf
)
:
	cyclicFvPatchField<Type>(ptf),
	subTimeLevel_(ptf.subTimeLevel_),
	phaseLagCyclicPatch_(ptf.phaseLagCyclicPatch_)
{}


template<class Type>
Foam::phaseLagCyclicFvPatchField<Type>::phaseLagCyclicFvPatchField
(
    const phaseLagCyclicFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
	cyclicFvPatchField<Type>(ptf, iF),
	subTimeLevel_(ptf.subTimeLevel_),
	phaseLagCyclicPatch_(ptf.phaseLagCyclicPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::phaseLagCyclicFvPatchField<Type>::phaseLaggedField() const
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

	tmp<Field<Type>> tPhaseLag(new Field<Type>(this->patchInternalField()));

	const dictionary& solDict = this->internalField().mesh().solutionDict();

	const label nT = solDict.subDict("harmonicBalance").getOrDefault<label>("instantsNumber",3);

	const RectangularMatrix<complex>& E = HB[HBZoneInstance].E();
	const RectangularMatrix<complex>& E_1 = HB[HBZoneInstance].EInv();

	const word& fieldName = this->internalField().name();

	const label& patchIndex = this->patch().index();
	const label& neighPatchIndex = this->phaseLagCyclicPatch().neighbPatchID();

	PtrList<Field<Type>> perioFields(nT);

	forAll(perioFields,i)
	{
		word itrName = Foam::name(i);
		word timeLevel = "subTimeLevel" + itrName;

		const objectRegistry& subLeveli = allSubLevels.lookupObject<objectRegistry>(timeLevel);

		if (subLeveli.found(fieldName))
		{
			const GeometricField<Type,fvPatchField,volMesh>& subTimeField =
					subLeveli.lookupObject<GeometricField<Type,fvPatchField,volMesh>>(fieldName);

			if (!subTimeField.boundaryField().operator()(patchIndex)
				|| !subTimeField.boundaryField().operator()(neighPatchIndex))
			{
				return tPhaseLag;
			}

		    const Field<Type>& iField = subTimeField.primitiveField();
		    const labelUList& nbrFaceCells =
		        this->phaseLagCyclicPatch().neighbFvPatch().faceCells();

		    Field<Type> subTimeInternali(this->size());

	        forAll(*this, facei)
	        {
	            subTimeInternali[facei] = iField[nbrFaceCells[facei]];
	        }

			perioFields.set
			(
				i,
				subTimeInternali
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

	Field<Type>& phaseLag = tPhaseLag.ref();

	phaseLag = Zero;

	forAll(phaseLag, facei)
	{
		forAll(perioFields, J)
		{
			phaseLag[facei] += D[subTimeLevel_][J]*perioFields[J][facei];
		}
	}

	return tPhaseLag;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::phaseLagCyclicFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->primitiveField();
    const labelUList& nbrFaceCells =
        this->phaseLagCyclicPatch().neighbFvPatch().faceCells();

    tmp<Field<Type>> tpnf(new Field<Type>(this->size()));
    Field<Type>& pnf = tpnf.ref();

    const word& fieldName = this->internalField().name();

    if (fieldName.startsWith("U") || fieldName == "p" || fieldName == "rho"
    		|| fieldName == "(he+(0.5*magSqr(U)))" || fieldName == "((he+(0.5*magSqr(U)))+(p|rho))"
    		|| fieldName == "sqrt((gamma|thermo:psi))" || fieldName == "c" || fieldName == "gamma"
    		|| fieldName == "H")
	{
	    Field<Type> jf(this->phaseLaggedField());

		forAll(*this, facei)
		{
			pnf[facei] = jf[facei];
		}
	}
	else
	{
		//Is it correct to transform the increment fields?
	    if (this->doTransform())
	    {
	        forAll(*this, facei)
	        {
	            pnf[facei] = transform
	            (
	                this->forwardT()[0], iField[nbrFaceCells[facei]]
	            );
	        }
	    }
	    else
	    {
	        forAll(*this, facei)
	        {
	            pnf[facei] = iField[nbrFaceCells[facei]];
	        }
	    }
	}

    return tpnf;
}


template<class Type>
void Foam::phaseLagCyclicFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    this->writeEntry("value", os);
}


// ************************************************************************* //
