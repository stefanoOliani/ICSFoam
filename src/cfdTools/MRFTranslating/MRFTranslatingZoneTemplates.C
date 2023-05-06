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

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"
#include "MRFTranslatingZone.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::MRFTranslatingZone::makeRelativeRhoFlux
(
    const RhoFieldType& rho,
    surfaceScalarField& phi
) const
{
    if (!active_)
    {
        return;
    }

    const surfaceVectorField& Sf = mesh_.Sf();

    const vector Velocity = frameVel_->value(mesh_.time().timeOutputValue());

    const vectorField& Sfi = Sf;
    scalarField& phii = phi.primitiveFieldRef();

    // Internal faces
    forAll(internalFaces_, i)
    {
        label facei = internalFaces_[i];
        phii[facei] -= rho[facei]*Velocity & Sfi[facei];
    }

    makeRelativeRhoFlux(rho.boundaryField(), phi.boundaryFieldRef());
}


template<class RhoFieldType>
void Foam::MRFTranslatingZone::makeRelativeRhoFlux
(
    const RhoFieldType& rho,
    FieldField<fvsPatchField, scalar>& phi
) const
{
    if (!active_)
    {
        return;
    }

    const surfaceVectorField& Sf = mesh_.Sf();

    const vector Velocity = frameVel_->value(mesh_.time().timeOutputValue());

    // Included patches
    forAll(includedFaces_, patchi)
    {
        forAll(includedFaces_[patchi], i)
        {
            label patchFacei = includedFaces_[patchi][i];

            phi[patchi][patchFacei] = 0.0;
        }
    }

    // Excluded patches
    forAll(excludedFaces_, patchi)
    {
        forAll(excludedFaces_[patchi], i)
        {
            label patchFacei = excludedFaces_[patchi][i];

            phi[patchi][patchFacei] -=
                rho[patchi][patchFacei] * Velocity & Sf.boundaryField()[patchi][patchFacei];
        }
    }
}


template<class RhoFieldType>
void Foam::MRFTranslatingZone::makeRelativeRhoFlux
(
    const RhoFieldType& rho,
    Field<scalar>& phi,
    const label patchi
) const
{
    if (!active_)
    {
        return;
    }

    const surfaceVectorField& Sf = mesh_.Sf();

    const vector Velocity = frameVel_->value(mesh_.time().timeOutputValue());

    // Included patches
    forAll(includedFaces_[patchi], i)
    {
        label patchFacei = includedFaces_[patchi][i];

        phi[patchFacei] = 0.0;
    }

    // Excluded patches
    forAll(excludedFaces_[patchi], i)
    {
        label patchFacei = excludedFaces_[patchi][i];

        phi[patchFacei] -=
            rho[patchFacei] * Velocity & Sf.boundaryField()[patchi][patchFacei];
    }
}


template<class RhoFieldType>
void Foam::MRFTranslatingZone::makeAbsoluteRhoFlux
(
    const RhoFieldType& rho,
    surfaceScalarField& phi
) const
{
    if (!active_)
    {
        return;
    }

    const surfaceVectorField& Sf = mesh_.Sf();

    const vector Velocity = frameVel_->value(mesh_.time().timeOutputValue());

    const vectorField& Sfi = Sf;
    scalarField& phii = phi.primitiveFieldRef();

    // Internal faces
    forAll(internalFaces_, i)
    {
        label facei = internalFaces_[i];
        phii[facei] += rho[facei]*Velocity & Sfi[facei];
    }

    surfaceScalarField::Boundary& phibf = phi.boundaryFieldRef();


    // Included patches
    forAll(includedFaces_, patchi)
    {
        forAll(includedFaces_[patchi], i)
        {
            label patchFacei = includedFaces_[patchi][i];

            phibf[patchi][patchFacei] +=
                rho.boundaryField()[patchi][patchFacei] * Velocity & Sf.boundaryField()[patchi][patchFacei];
        }
    }

    // Excluded patches
    forAll(excludedFaces_, patchi)
    {
        forAll(excludedFaces_[patchi], i)
        {
            label patchFacei = excludedFaces_[patchi][i];

            phibf[patchi][patchFacei] +=
                rho.boundaryField()[patchi][patchFacei] * Velocity & Sf.boundaryField()[patchi][patchFacei];
        }
    }
}


template<class Type>
void Foam::MRFTranslatingZone::zero
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& phi
) const
{
    if (!active_)
    {
        return;
    }

    Field<Type>& phii = phi.primitiveFieldRef();

    forAll(internalFaces_, i)
    {
        phii[internalFaces_[i]] = Zero;
    }

    typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary& phibf =
        phi.boundaryFieldRef();

    forAll(includedFaces_, patchi)
    {
        forAll(includedFaces_[patchi], i)
        {
            phibf[patchi][includedFaces_[patchi][i]] = Zero;
        }
    }

    forAll(excludedFaces_, patchi)
    {
        forAll(excludedFaces_[patchi], i)
        {
            phibf[patchi][excludedFaces_[patchi][i]] = Zero;
        }
    }
}


// ************************************************************************* //
