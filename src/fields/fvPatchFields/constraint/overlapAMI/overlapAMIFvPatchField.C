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

#include "overlapAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::overlapAMIFvPatchField<Type>::overlapAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    overlapAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF),
    overlapAMIPatch_(refCast<const overlapAMIFvPatch>(p))
{}


template<class Type>
Foam::overlapAMIFvPatchField<Type>::overlapAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    overlapAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF, dict, dict.found("value")),
    overlapAMIPatch_(refCast<const overlapAMIFvPatch>(p, dict))
{
    if (!isA<overlapAMIFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    if (!dict.found("value") && this->coupled())
    {
        this->evaluate(Pstream::commsTypes::blocking);
    }
}


template<class Type>
Foam::overlapAMIFvPatchField<Type>::overlapAMIFvPatchField
(
    const overlapAMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    overlapAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    overlapAMIPatch_(refCast<const overlapAMIFvPatch>(p))
{
    if (!isA<overlapAMIFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::overlapAMIFvPatchField<Type>::overlapAMIFvPatchField
(
    const overlapAMIFvPatchField<Type>& ptf
)
:
    overlapAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    overlapAMIPatch_(ptf.overlapAMIPatch_)
{}


template<class Type>
Foam::overlapAMIFvPatchField<Type>::overlapAMIFvPatchField
(
    const overlapAMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    overlapAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    overlapAMIPatch_(ptf.overlapAMIPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::overlapAMIFvPatchField<Type>::coupled() const
{
    return overlapAMIPatch_.coupled();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::overlapAMIFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->primitiveField();
    const labelUList& nbrFaceCells =
        overlapAMIPatch_.overlapAMIPatch().neighbPatch().faceCells();

    Field<Type> pnf(iField, nbrFaceCells);

    const word& fieldName = this->internalField().name();

    tmp<Field<Type>> tpnf;
    if (overlapAMIPatch_.applyLowWeightCorrection())
    {
    		tpnf = overlapAMIPatch_.interpolate(pnf, this->patchInternalField()());
    }
    else
    {
    	if (fieldName.startsWith("U"))
		{
    		if (fieldName == "U")
    		{
    			tpnf = overlapAMIPatch_.interpolate(pnf);
    		}
    		else
    		{
    			if (fieldName == "U.component(0)")
    			{
    				tpnf = overlapAMIPatch_.interpolate(pnf, 0);
    			}
    			else if (fieldName == "U.component(1)")
				{
    				tpnf = overlapAMIPatch_.interpolate(pnf, 1);
				}
    			else if  (fieldName == "U.component(2)")
    			{
    				tpnf = overlapAMIPatch_.interpolate(pnf, 2);
    			}
    			else
    			{
    				tpnf = overlapAMIPatch_.untransfInterp(pnf);
    			}
    		}
		}
		else
		{
			tpnf = overlapAMIPatch_.interpolate(pnf);
		}
    }

    return tpnf;
}


template<class Type>
const Foam::overlapAMIFvPatchField<Type>&
Foam::overlapAMIFvPatchField<Type>::neighbourPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->primitiveField()
        );

    return refCast<const overlapAMIFvPatchField<Type>>
    (
        fld.boundaryField()[overlapAMIPatch_.neighbPatchID()]
    );
}


template<class Type>
void Foam::overlapAMIFvPatchField<Type>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
	NotImplemented;
}


template<class Type>
void Foam::overlapAMIFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{

    const labelList& nbrFaceCells =
                    lduAddr.patchAddr
                    (
                        overlapAMIPatch_.overlapAMIPatch().neighbPatchID()
                    );

    Field<Type> pnf(psiInternal, nbrFaceCells);

    if (overlapAMIPatch_.applyLowWeightCorrection())
    {
        Field<Type> pif(psiInternal, overlapAMIPatch_.faceCells());
        pnf = overlapAMIPatch_.interpolate(pnf, pif);
    }
    else
    {
        pnf = overlapAMIPatch_.interpolate(pnf);
    }

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}


template<class Type>
void Foam::overlapAMIFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
