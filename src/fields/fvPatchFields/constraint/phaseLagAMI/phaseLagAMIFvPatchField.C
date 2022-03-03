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

#include "phaseLagAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::phaseLagAMIFvPatchField<Type>::phaseLagAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    phaseLagAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF),
    phaseLagAMIPatch_(refCast<const phaseLagAMIFvPatch>(p))
{}


template<class Type>
Foam::phaseLagAMIFvPatchField<Type>::phaseLagAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    phaseLagAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF, dict, dict.found("value")),
    phaseLagAMIPatch_(refCast<const phaseLagAMIFvPatch>(p, dict))
{
    if (!isA<phaseLagAMIFvPatch>(p))
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
Foam::phaseLagAMIFvPatchField<Type>::phaseLagAMIFvPatchField
(
    const phaseLagAMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    phaseLagAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    phaseLagAMIPatch_(refCast<const phaseLagAMIFvPatch>(p))
{
    if (!isA<phaseLagAMIFvPatch>(this->patch()))
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
Foam::phaseLagAMIFvPatchField<Type>::phaseLagAMIFvPatchField
(
    const phaseLagAMIFvPatchField<Type>& ptf
)
:
    phaseLagAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    phaseLagAMIPatch_(ptf.phaseLagAMIPatch_)
{}


template<class Type>
Foam::phaseLagAMIFvPatchField<Type>::phaseLagAMIFvPatchField
(
    const phaseLagAMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    phaseLagAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    phaseLagAMIPatch_(ptf.phaseLagAMIPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::phaseLagAMIFvPatchField<Type>::coupled() const
{
    return phaseLagAMIPatch_.coupled();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::phaseLagAMIFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->primitiveField();
    const labelUList& nbrFaceCells =
        phaseLagAMIPatch_.phaseLagAMIPatch().neighbPatch().faceCells();

    Field<Type> pnf(iField, nbrFaceCells);

    const word& fieldName = this->internalField().name();

    tmp<Field<Type>> tpnf;
    if (phaseLagAMIPatch_.applyLowWeightCorrection())
    {
    	if (fieldName.startsWith("U") || fieldName == "p" || fieldName == "rho"
			|| fieldName == "(he+(0.5*magSqr(U)))" || fieldName == "((he+(0.5*magSqr(U)))+(p|rho))"
			|| fieldName == "sqrt((gamma|thermo:psi))" || fieldName == "c" || fieldName == "gamma"
			|| fieldName == "H")
    	{
    		tpnf = phaseLagAMIPatch_.interpolate(pnf, fieldName, this->patchInternalField()());
    	}
    	else
    	{
            tpnf = phaseLagAMIPatch_.interpolate(pnf, this->patchInternalField()());
    	}
    }
    else
    {
    	if (fieldName.startsWith("U") || fieldName == "p" || fieldName == "rho"
			|| fieldName == "(he+(0.5*magSqr(U)))" || fieldName == "((he+(0.5*magSqr(U)))+(p|rho))"
			|| fieldName == "sqrt((gamma|thermo:psi))" || fieldName == "c" || fieldName == "gamma"
			|| fieldName == "H")
    	{
    		tpnf = phaseLagAMIPatch_.interpolate(pnf, fieldName);
    	}
    	else
    	{
    		tpnf = phaseLagAMIPatch_.interpolate(pnf);
    	}
    }

    //Transform should have been already performed during interpolate
//    if (doTransform())
//    {
//        tpnf.ref() = transform(forwardT(), tpnf());
//    }

    return tpnf;
}


template<class Type>
const Foam::phaseLagAMIFvPatchField<Type>&
Foam::phaseLagAMIFvPatchField<Type>::neighbourPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->primitiveField()
        );

    return refCast<const phaseLagAMIFvPatchField<Type>>
    (
        fld.boundaryField()[phaseLagAMIPatch_.neighbPatchID()]
    );
}


template<class Type>
void Foam::phaseLagAMIFvPatchField<Type>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
	NotImplemented;
}


template<class Type>
void Foam::phaseLagAMIFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{

	Info<<"Should not call it"<<endl;
    const labelUList& nbrFaceCells =
        phaseLagAMIPatch_.phaseLagAMIPatch().neighbPatch().faceCells();

    Field<Type> pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf);

    const word& fieldName = this->internalField().name();

    if (phaseLagAMIPatch_.applyLowWeightCorrection())
    {
        Field<Type> pif(psiInternal, phaseLagAMIPatch_.faceCells());
        pnf = phaseLagAMIPatch_.interpolate(pnf, fieldName, pif);
    }
    else
    {
        pnf = phaseLagAMIPatch_.interpolate(pnf, fieldName);
    }

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, coeffs, pnf);
}


template<class Type>
void Foam::phaseLagAMIFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
