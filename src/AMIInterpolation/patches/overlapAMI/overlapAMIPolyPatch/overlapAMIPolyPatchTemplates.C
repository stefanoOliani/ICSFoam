/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

    Hrvoje Jasak, Wikki Ltd.  All rights reserved
    Fethi Tekin, All rights reserved.
    Oliver Borm, All rights reserved.

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::overlapAMIPolyPatch::expandData(const Field<Type>& pf) const
{
    // Check and expand the field from patch size to zone size
    if (pf.size() != this->size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > overlapAMIPolyPatch::expandData"
            "("
            "    const Field<Type>& pf"
            ") const"
        )   << "Incorrect patch field size.  Field size: "
            << pf.size() << " patch size: " << this->size()
            << abort(FatalError);
    }

    const label ncp = nCopies();

	const scalar myAngle = 360.0/scalar(ncp);

    tmp<Field<Type> > texpandField
    (
	   new Field<Type>(ncp*pf.size(), pTraits<Type>::zero)
    );

    Field<Type>& expandField = texpandField.ref();

    for (label copyI = 0; copyI < ncp; copyI++)
    {
    	// Calculate transform
		const tensor curRotation = this->RodriguesRotation(rotationAxis_, copyI*myAngle);

		const label offset = copyI*pf.size();

		forAll (pf, faceI)
		{
			 const label zId = this->whichFace(this->start() + faceI);
			 expandField[offset + zId] = Foam::transform(curRotation, pf[faceI]);
		}
    }

    return texpandField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::overlapAMIPolyPatch::expandData(const Field<Type>& pf, label cmpt) const
{
	Info<<"should not call it"<<endl;

    // Check and expand the field from patch size to zone size
    if (pf.size() != this->size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > overlapAMIPolyPatch::expandData"
            "("
            "    const Field<Type>& pf"
            ") const"
        )   << "Incorrect patch field size.  Field size: "
            << pf.size() << " patch size: " << this->size()
            << abort(FatalError);
    }

    const label ncp = nCopies();

	const scalar myAngle = 360.0/scalar(ncp);

    tmp<Field<Type> > texpandField
    (
	   new Field<Type>(ncp*pf.size(), pTraits<Type>::zero)
    );

    Field<Type>& expandField = texpandField.ref();

    for (label copyI = 0; copyI < ncp; copyI++)
    {
    	// Calculate transform
		const tensor curRotation = this->RodriguesRotation(rotationAxis_, copyI*myAngle);

		const label offset = copyI*pf.size();

		forAll (pf, faceI)
		{
			 const label zId = this->whichFace(this->start() + faceI);
			 expandField[offset + zId] = Foam::transform(curRotation, pf[faceI]);
		}
    }

    return texpandField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::overlapAMIPolyPatch::untransfExpandData(const Field<Type>& pf) const
{
    // Check and expand the field from patch size to zone size
    if (pf.size() != this->size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > overlapAMIPolyPatch::expandData"
            "("
            "    const Field<Type>& pf"
            ") const"
        )   << "Incorrect patch field size.  Field size: "
            << pf.size() << " patch size: " << this->size()
            << abort(FatalError);
    }

    const label ncp = nCopies();

    tmp<Field<Type> > texpandField
    (
	   new Field<Type>(ncp*pf.size(), pTraits<Type>::zero)
    );

    Field<Type>& expandField = texpandField.ref();

    for (label copyI = 0; copyI < ncp; copyI++)
    {
		const label offset = copyI*pf.size();

		forAll (pf, faceI)
		{
			 const label zId = this->whichFace(this->start() + faceI);
			 expandField[offset + zId] =pf[faceI];
		}
    }

    return texpandField;
}


template<class Type>
Foam::UList<Type>
Foam::overlapAMIPolyPatch::expandData(const UList<Type>& defaultValues) const
{
	if (defaultValues.size())
	{
		// Check and expand the field from patch size to zone size
		if (defaultValues.size() != this->size())
		{
			FatalErrorIn
			(
				"tmp<Field<Type> > overlapAMIPolyPatch::expandData"
				"("
				"    UList<Type>& defaultValues"
				") const"
			)   << "Incorrect patch field size.  Field size: "
				<< defaultValues.size() << " patch size: " << this->size()
				<< abort(FatalError);
		}

		const label ncp = nCopies();

		const scalar myAngle = 360.0/scalar(ncp);

		Type dfl = *defaultValues.cdata();
		Type* dflPtr = &dfl;

		UList<Type> expandField(dflPtr, nCopies_*defaultValues.size());

		for (label copyI = 0; copyI < ncp; copyI++)
		{
			// Calculate transform
			if (owner())
			{
				 const tensor curRotation = this->RodriguesRotation(rotationAxis_, copyI*myAngle);

				 const label offset = copyI*defaultValues.size();

				 forAll (defaultValues, faceI)
				 {
					 const label zId = this->whichFace(this->start() + faceI);
					 expandField[offset + zId] = Foam::transform(curRotation, defaultValues[faceI]);
				 }
			}
			else
			{
				const tensor curRotation = this->RodriguesRotation(rotationAxis_, copyI*myAngle);

				const label offset = copyI*defaultValues.size();

				 forAll (defaultValues, faceI)
				 {
					 const label zId = this->whichFace(this->start() + faceI);
					 expandField[offset + zId] = Foam::transform(curRotation, defaultValues[faceI]);
				 }
			}
		}

		return expandField;
	}

	return defaultValues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::overlapAMIPolyPatch::interpolate
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
Foam::tmp<Foam::Field<Type>> Foam::overlapAMIPolyPatch::interpolate
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolate(tFld(), defaultValues);
}


template<class Type, class CombineOp>
void Foam::overlapAMIPolyPatch::interpolate
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


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::overlapAMIPolyPatch::interpolate
(
    const Field<Type>& fld,
	label cmpt,
    const UList<Type>& defaultValues
) const
{
	// Expand data
	tmp<Field<Type> > expandDataTmp = neighbPatch().expandData(fld, cmpt);
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
Foam::tmp<Foam::Field<Type>> Foam::overlapAMIPolyPatch::interpolate
(
    const tmp<Field<Type>>& tFld,
	label cmpt,
    const UList<Type>& defaultValues
) const
{
    return interpolate(tFld(), cmpt, defaultValues);
}


template<class Type, class CombineOp>
void Foam::overlapAMIPolyPatch::interpolate
(
    const UList<Type>& fld,
	label cmpt,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
	// Expand data
	tmp<Field<Type> > expandDataTmp = neighbPatch().expandData(fld, cmpt);

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
Foam::tmp<Foam::Field<Type>> Foam::overlapAMIPolyPatch::untransfInterp
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
	// Expand data
	tmp<Field<Type> > expandDataTmp = neighbPatch().untransfExpandData(fld);
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
Foam::tmp<Foam::Field<Type>> Foam::overlapAMIPolyPatch::untransfInterp
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return untransfInterp(tFld(), defaultValues);
}


template<class Type, class CombineOp>
void Foam::overlapAMIPolyPatch::untransfInterp
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
	// Expand data
	tmp<Field<Type> > expandDataTmp = neighbPatch().untransfExpandData(fld);

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
