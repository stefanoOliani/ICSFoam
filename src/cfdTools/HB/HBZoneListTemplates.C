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

#include "HBZoneList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::HBZoneList::addingSource
(
	PtrList<Field<Type>>& source,
	PtrList<PtrList<GeometricField<Type, fvPatchField, volMesh>>>& vars,
	label varNum
) const
{
    forAll(*this, i)
    {
        operator[](i).addSource(source, vars, varNum);
    }
}


template<class Type>
void Foam::HBZoneList::addingSource
(
	fvMatrix<Type>& eqn,
	PtrList<GeometricField<Type, fvPatchField, volMesh>>& fieldPtr,
	label instantNo
) const
{
    forAll(*this, i)
    {
        operator[](i).addSource(eqn, fieldPtr, instantNo);
    }
}

template<class Type>
void Foam::HBZoneList::addingSource
(
	const volScalarField& rho,
	fvMatrix<Type>& eqn,
	PtrList<GeometricField<Type, fvPatchField, volMesh>>& fieldPtr,
	label instantNo
) const
{
    forAll(*this, i)
    {
        operator[](i).addSource(rho, eqn, fieldPtr, instantNo);
    }
}


template<class Type>
void Foam::HBZoneList::factStep
(
	PtrList<GeometricField<Type, fvPatchField, volMesh>>& fieldPtr,
	PtrList<scalarField>& deltaTField
) const
{
    forAll(*this, i)
    {
        operator[](i).factorizationStep(fieldPtr, deltaTField);
    }
}


template<class Type>
void Foam::HBZoneList::reconstruct
(
	GeometricField<Type, fvPatchField, volMesh>& reconstrFld,
	PtrList<GeometricField<Type, fvPatchField, volMesh>>& fieldPtr,
	scalar actTime
) const
{
    forAll(*this, i)
    {
        operator[](i).reconstruct(reconstrFld, fieldPtr, actTime);
    }
}


// ************************************************************************* //
