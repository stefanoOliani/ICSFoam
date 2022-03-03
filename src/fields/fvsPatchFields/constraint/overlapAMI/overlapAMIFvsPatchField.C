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

#include "overlapAMIFvsPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::overlapAMIFvsPatchField<Type>::overlapAMIFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(p, iF),
    overlapAMIPatch_(refCast<const overlapAMIFvPatch>(p))
{}


template<class Type>
Foam::overlapAMIFvsPatchField<Type>::overlapAMIFvsPatchField
(
    const overlapAMIFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvsPatchField<Type>(ptf, p, iF, mapper),
    overlapAMIPatch_(refCast<const overlapAMIFvPatch>(p))
{
    if (!isA<overlapAMIFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::overlapAMIFvsPatchField<Type>::overlapAMIFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    coupledFvsPatchField<Type>(p, iF, dict),
    overlapAMIPatch_(refCast<const overlapAMIFvPatch>(p, dict))
{
    if (!isA<overlapAMIFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "patch " << this->patch().index() << " not overlapAMI type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::overlapAMIFvsPatchField<Type>::overlapAMIFvsPatchField
(
    const overlapAMIFvsPatchField<Type>& ptf
)
:
    coupledFvsPatchField<Type>(ptf),
    overlapAMIPatch_(ptf.overlapAMIPatch_)
{}


template<class Type>
Foam::overlapAMIFvsPatchField<Type>::overlapAMIFvsPatchField
(
    const overlapAMIFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(ptf, iF),
    overlapAMIPatch_(ptf.overlapAMIPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::overlapAMIFvsPatchField<Type>::coupled() const
{
    if
    (
        Pstream::parRun()
     || (
            this->overlapAMIPatch_.size()
         && this->overlapAMIPatch_.overlapAMIPatch().neighbPatch().size()
        )
    )
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
