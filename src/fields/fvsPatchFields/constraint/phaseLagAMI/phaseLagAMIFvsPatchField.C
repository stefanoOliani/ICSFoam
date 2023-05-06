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

#include "phaseLagAMIFvsPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::phaseLagAMIFvsPatchField<Type>::phaseLagAMIFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(p, iF),
    phaseLagAMIPatch_(refCast<const phaseLagAMIFvPatch>(p))
{}


template<class Type>
Foam::phaseLagAMIFvsPatchField<Type>::phaseLagAMIFvsPatchField
(
    const phaseLagAMIFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvsPatchField<Type>(ptf, p, iF, mapper),
    phaseLagAMIPatch_(refCast<const phaseLagAMIFvPatch>(p))
{
    if (!isA<phaseLagAMIFvPatch>(this->patch()))
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
Foam::phaseLagAMIFvsPatchField<Type>::phaseLagAMIFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    coupledFvsPatchField<Type>(p, iF, dict),
    phaseLagAMIPatch_(refCast<const phaseLagAMIFvPatch>(p, dict))
{
    if (!isA<phaseLagAMIFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "patch " << this->patch().index() << " not phaseLagAMI type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::phaseLagAMIFvsPatchField<Type>::phaseLagAMIFvsPatchField
(
    const phaseLagAMIFvsPatchField<Type>& ptf
)
:
    coupledFvsPatchField<Type>(ptf),
    phaseLagAMIPatch_(ptf.phaseLagAMIPatch_)
{}


template<class Type>
Foam::phaseLagAMIFvsPatchField<Type>::phaseLagAMIFvsPatchField
(
    const phaseLagAMIFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(ptf, iF),
    phaseLagAMIPatch_(ptf.phaseLagAMIPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::phaseLagAMIFvsPatchField<Type>::coupled() const
{
    if
    (
        Pstream::parRun()
     || (
            this->phaseLagAMIPatch_.size()
         && this->phaseLagAMIPatch_.phaseLagAMIPatch().neighbPatch().size()
        )
    )
    {
        return true;
    }

    return false;
}


// ************************************************************************* //
