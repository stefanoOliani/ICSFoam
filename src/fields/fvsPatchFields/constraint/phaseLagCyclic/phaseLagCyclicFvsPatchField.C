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

#include "phaseLagCyclicFvsPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::phaseLagCyclicFvsPatchField<Type>::phaseLagCyclicFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    cyclicFvsPatchField<Type>(p, iF),
    phaseLagCyclicPatch_(refCast<const phaseLagCyclicFvPatch>(p))
{}


template<class Type>
Foam::phaseLagCyclicFvsPatchField<Type>::phaseLagCyclicFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    cyclicFvsPatchField<Type>(p, iF, dict),
    phaseLagCyclicPatch_(refCast<const phaseLagCyclicFvPatch>(p, dict))
{
    if (!isA<phaseLagCyclicFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "patch " << this->patch().index() << " not phaseLagCyclic type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::phaseLagCyclicFvsPatchField<Type>::phaseLagCyclicFvsPatchField
(
    const phaseLagCyclicFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicFvsPatchField<Type>(ptf, p, iF, mapper),
    phaseLagCyclicPatch_(refCast<const phaseLagCyclicFvPatch>(p))
{
    if (!isA<phaseLagCyclicFvPatch>(this->patch()))
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
Foam::phaseLagCyclicFvsPatchField<Type>::phaseLagCyclicFvsPatchField
(
    const phaseLagCyclicFvsPatchField<Type>& ptf
)
:
    cyclicFvsPatchField<Type>(ptf),
    phaseLagCyclicPatch_(ptf.phaseLagCyclicPatch_)
{}


template<class Type>
Foam::phaseLagCyclicFvsPatchField<Type>::phaseLagCyclicFvsPatchField
(
    const phaseLagCyclicFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    cyclicFvsPatchField<Type>(ptf, iF),
    phaseLagCyclicPatch_(ptf.phaseLagCyclicPatch_)
{}


// ************************************************************************* //
