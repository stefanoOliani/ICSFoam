/*---------------------------------------------------------------------------*\

    ICSFoam: a library for Implicit Coupled Simulations in OpenFOAM
  
    Copyright (C) 2022	Stefano Oliani

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

#ifndef cellFaceFunctions_H
#define cellFaceFunctions_H

#include "geometricOneField.H"

// Select from one or the other surface field comparing to a threshold value
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
surfaceFieldSelect
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfL,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfR,
    const surfaceScalarField& select,
    const scalar thres
)
{
    const fvMesh& mesh = sfL.mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "select("+sfL.name()+","+sfR.name()+')',
                sfL.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            sfL.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf.ref();

    if (sfL.dimensions() != sfR.dimensions())
    {
        FatalErrorInFunction << "Dimensions of sfL and sfR differ: " << nl
                     << sfL.dimensions() << " != " << sfL.dimensions() << nl
                     << exit(FatalError);
    }

    forAll(select, facei)
    {
        sf[facei] = (select[facei] >= thres ? sfL[facei] : sfR[facei]);
    }

    forAll(select.boundaryField(), patchi)
    {
        Field<Type>& psf = sf.boundaryFieldRef()[patchi];
        const Field<Type>& psfL = sfL.boundaryField()[patchi];
        const Field<Type>& psfR = sfR.boundaryField()[patchi];
        const scalarField& pSel = select.boundaryField()[patchi];
        forAll(pSel, facei)
        {
            psf[facei] = (pSel[facei] >= thres ? psfL[facei] : psfR[facei]);   
        }
    }

    return tsf;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
surfaceFieldSelect
(
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& sfL,
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& sfR,
    const surfaceScalarField& select,
    const scalar thres
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf =
        surfaceFieldSelect(sfL(), sfR(), select, thres);
    sfL.clear();
    sfR.clear();
    return tsf;
}


// Call a function for each face of a field; return result as a field
template<class Type, class Lambda>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
calcForEachFace
(
    Lambda fn,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
)
{
    const fvMesh& mesh = sf.mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfout
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "calc("+sf.name()+')',
                sf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            sf.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sfout = tsfout.ref();

    forAll(sfout, facei)
    {
        sfout[facei] = fn(sf[facei]);
    }

    forAll(sfout.boundaryField(), patchi)
    {
        Field<Type>& psfout = sfout.boundaryFieldRef()[patchi];
        const Field<Type>& psf = sf.boundaryField()[patchi];
        forAll(psfout, facei)
        {
            psfout[facei] = fn(psf[facei]);
        }
    }

    return tsfout;
}


// Call a function for each face of a field; return result as a field
template<class Type, class Lambda>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
calcForEachFace
(
    Lambda fn,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf1,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf2   
)
{
    const fvMesh& mesh = sf1.mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfout
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "calc("+sf1.name()+','+sf2.name()+')',
                sf1.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            sf1.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sfout = tsfout.ref();

    forAll(sfout, facei)
    {
        sfout[facei] = fn(sf1[facei], sf2[facei]);
    }

    forAll(sfout.boundaryField(), patchi)
    {
        Field<Type>& psfout = sfout.boundaryFieldRef()[patchi];
        const Field<Type>& psf1 = sf1.boundaryField()[patchi];
        const Field<Type>& psf2 = sf2.boundaryField()[patchi];
        forAll(psfout, facei)
        {
            psfout[facei] = fn(psf1[facei], psf2[facei]);
        }
    }

    return tsfout;
}

#endif

// ************************************************************************* //

