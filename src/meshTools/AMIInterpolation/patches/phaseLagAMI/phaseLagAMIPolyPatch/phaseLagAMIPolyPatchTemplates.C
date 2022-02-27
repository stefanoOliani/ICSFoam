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

#include "IOHBZoneList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::phaseLagAMIPolyPatch::expandData(const Field<Type>& pf) const
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

    // Get the periodic patch
    const coupledPolyPatch& periodicPatch
    (
        refCast<const coupledPolyPatch>
        (
            boundaryMesh()[periodicPatchID()]
        )
    );

    tmp<Field<Type>> texpandField
    (
        new Field<Type>(pf)
    );

    Field<Type>& expandField = texpandField.ref();
	Field<Type> transfField(pf);

    for (label copyI = 0; copyI < nTransformsBwd_; copyI++)
    {
        // Calculate transform
		 const tensorField& curTransform = periodicPatch.reverseT();

    	if (curTransform.size())
    	{
    		 transfField = Foam::transform(curTransform, transfField);
    	}

		expandField.append(transfField);
    }

    transfField = pf;

    for (label copyI = 0; copyI < nTransformsFwd_; copyI++)
    {
        // Calculate transform
		 const tensorField& curTransform = periodicPatch.forwardT();

    	if (curTransform.size())
    	{
    		 transfField = Foam::transform(curTransform, transfField);
    	}

		 expandField.append(transfField);
    }

    return texpandField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::phaseLagAMIPolyPatch::expandData(const Field<Type>& pf, const word& fieldName) const
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

	PtrList<Field<Type>> perioFields(nT);

    tmp<Field<Type>> texpandField
    (
        new Field<Type>(pf)
    );

	forAll(perioFields,i)
	{
		word itrName = Foam::name(i);
		word timeLevel = "subTimeLevel" + itrName;

		const objectRegistry& subLeveli = allSubLevels.lookupObject<objectRegistry>(timeLevel);

		if(subLeveli.found(fieldName))
		{
			const GeometricField<Type,fvPatchField,volMesh>& subTimeField =
					subLeveli.lookupObject<GeometricField<Type,fvPatchField,volMesh>>(fieldName);

			if (!subTimeField.boundaryField().operator()(patchIndex)
				|| !subTimeField.boundaryField().operator()(neighPatchIndex))
			{
				const Field<Type> subTimeInternali = pf;

				perioFields.set
				(
					i,
					subTimeInternali
				);
			}
			else
			{
				const Field<Type> subTimeInternali
										= subTimeField.boundaryField()[patchIndex].patchInternalField();

				perioFields.set
				(
					i,
					subTimeInternali
				);
			}
		}
		else
		{
			const Field<Type> subTimeInternali = pf;

			perioFields.set
			(
				i,
				subTimeInternali
			);
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

    Field<Type>& expandField = texpandField.ref();
	Field<Type> transfField(pf.size(), Zero);

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

		forAll(pf, facei)
		{
			forAll(perioFields, J)
			{
				transfField[facei] += D[subTimeLevel_][J]*perioFields[J][facei];
			}
		}

        // Calculate transform
		const tensorField& curTransform = periodicPatch.reverseT();

    	if (curTransform.size())
    	{
    		 transfField = Foam::transform(curTransform, transfField);
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

		forAll(pf, facei)
		{
			forAll(perioFields, J)
			{
				transfField[facei] += D[subTimeLevel_][J]*perioFields[J][facei];
			}
		}

        // Calculate transform
		const tensorField& curTransform = periodicPatch.forwardT();

    	if (curTransform.size())
    	{
    		 transfField = Foam::transform(curTransform, transfField);
    	}

		 expandField.append(transfField);

		 transfField = Zero;
    }

    return texpandField;
}


template<class Type>
Foam::UList<Type>
Foam::phaseLagAMIPolyPatch::expandData(const UList<Type>& defaultValues) const
{
	if (defaultValues.size())
	{
		// Check and expand the field from patch size to zone size
		if (defaultValues.size() != this->size())
		{
			FatalErrorIn
			(
				"tmp<Field<Type> > phaseLagAMIPolyPatch::expandData"
				"("
				"    UList<Type>& defaultValues"
				") const"
			)   << "Incorrect patch field size.  Field size: "
				<< defaultValues.size() << " patch size: " << this->size()
				<< abort(FatalError);
		}

		// Get the periodic patch
		const coupledPolyPatch& periodicPatch
		(
			refCast<const coupledPolyPatch>
			(
				boundaryMesh()[periodicPatchID()]
			)
		);


		Type dfl = *defaultValues.cdata();
		Type* dflPtr = &dfl;

		UList<Type> expandField(dflPtr, nTransforms_*defaultValues.size());

		for (label copyI = 0; copyI < nTransforms_; copyI++)
		{
			// Calculate transform
			if (owner())
			{
				 const tensorField& curTransform = periodicPatch.forwardT();

				 const label offset = copyI*defaultValues.size();

				 forAll (defaultValues, faceI)
				 {
					 const label zId = this->whichFace(this->start() + faceI);
					 expandField[offset + zId] = Foam::transform(curTransform[faceI], defaultValues[faceI]);
				 }
			}
			else
			{
				const tensorField& curTransform = periodicPatch.reverseT();

				const label offset = copyI*defaultValues.size();

				 forAll (defaultValues, faceI)
				 {
					 const label zId = this->whichFace(this->start() + faceI);
					 expandField[offset + zId] = Foam::transform(curTransform[faceI], defaultValues[faceI]);
				 }
			}
		}

		return expandField;
	}

	return defaultValues;
}



template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::phaseLagAMIPolyPatch::interpolate
(
    const Field<Type>& fld,
	const word& fieldName,
    const UList<Type>& defaultValues
) const
{
	// Expand data

	tmp<Field<Type> > expandDataTmp = neighbPatch().expandData(fld, fieldName);
	Field<Type>& expandData = expandDataTmp.ref();

	UList<Type>  expandDefault = neighbPatch().expandData(defaultValues);

	tmp<Field<Type> > tresult(new Field<Type>());
	Field<Type>& result = tresult.ref();

    if (owner())
    {
        result = AMI().interpolateToSource(expandData, expandDefault);
        // Truncate to size
        result.setSize(this->size());

        return tresult;
    }
    else
    {
        result = neighbPatch().AMI().interpolateToTarget(expandData, expandDefault);
        // Truncate to size
        result.setSize(this->size());

        return tresult;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::phaseLagAMIPolyPatch::interpolate
(
    const tmp<Field<Type>>& tFld,
	const word& fieldName,
    const UList<Type>& defaultValues
) const
{
    return interpolate(tFld(), fieldName, defaultValues);
}


template<class Type, class CombineOp>
void Foam::phaseLagAMIPolyPatch::interpolate
(
    const UList<Type>& fld,
	const word& fieldName,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
	// Expand data
	tmp<Field<Type> > expandDataTmp = neighbPatch().expandData(fld, fieldName);
	Field<Type>& expandData = expandDataTmp.ref();

	UList<Type> expandDefault = neighbPatch().expandData(defaultValues);

    if (owner())
    {
        AMI().interpolateToSource
        (
            expandData,
            cop,
            result,
	    expandDefault
        );

        // Truncate to size
	result.setSize(this->size());
    }
    else
    {
        neighbPatch().AMI().interpolateToTarget
        (
            expandData,
            cop,
            result,
	    expandDefault
        );

        // Truncate to size
	result.setSize(this->size());
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::phaseLagAMIPolyPatch::interpolate
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
	// Expand data
	tmp<Field<Type> > expandDataTmp = neighbPatch().expandData(fld);
	Field<Type>& expandData = expandDataTmp.ref();

	UList<Type>  expandDefault = neighbPatch().expandData(defaultValues);

	tmp<Field<Type> > tresult(new Field<Type>());
	Field<Type>& result = tresult.ref();

    if (owner())
    {
        result = AMI().interpolateToSource(expandData, expandDefault);
        // Truncate to size
        result.setSize(this->size());

        return tresult;
    }
    else
    {
        result = neighbPatch().AMI().interpolateToTarget(expandData, expandDefault);
        // Truncate to size
        result.setSize(this->size());

        return tresult;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::phaseLagAMIPolyPatch::interpolate
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolate(tFld(), defaultValues);
}


template<class Type, class CombineOp>
void Foam::phaseLagAMIPolyPatch::interpolate
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
	// Expand data
	tmp<Field<Type> > expandDataTmp = neighbPatch().expandData(fld);

	Field<Type>& expandData = expandDataTmp.ref();

	UList<Type> expandDefault = neighbPatch().expandData(defaultValues);

    if (owner())
    {
        AMI().interpolateToSource
        (
            expandData,
            cop,
            result,
	    expandDefault
        );

        // Truncate to size
	result.setSize(this->size());
    }
    else
    {
        neighbPatch().AMI().interpolateToTarget
        (
            expandData,
            cop,
            result,
	    expandDefault
        );

        // Truncate to size
	result.setSize(this->size());
    }
}


// ************************************************************************* //
