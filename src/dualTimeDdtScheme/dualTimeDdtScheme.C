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

#include "dualTimeDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
const volScalarField& dualTimeDdtScheme<Type>::localRPseudoDeltaT() const
{
    return mesh().objectRegistry::template lookupObject<volScalarField>
        (rPseudoDeltaTName_);
}

// No pseudotime contrib to fvc scheme, only fvm

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
dualTimeDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    return scheme_->fvcDdt(dt);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
dualTimeDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_->fvcDdt(vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
dualTimeDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_->fvcDdt(rho, vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
dualTimeDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_->fvcDdt(rho, vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
dualTimeDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_->fvcDdt(alpha, rho, vf);
}


template<class Type>
tmp<fvMatrix<Type> >
dualTimeDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm = scheme_->fvmDdt(vf);
    fvMatrix<Type>& fvm = tfvm.ref();
    const scalarField& rPseudoDeltaT = localRPseudoDeltaT().internalField();

    fvm.diag() += rPseudoDeltaT*mesh().Vsc();
    fvm.source() += rPseudoDeltaT*vf.internalField()*mesh().Vsc();

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
dualTimeDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm = scheme_->fvmDdt(rho,vf);
    fvMatrix<Type>& fvm = tfvm.ref();
    const scalarField& rPseudoDeltaT = localRPseudoDeltaT().internalField();

    fvm.diag() += rPseudoDeltaT*rho.value()*mesh().Vsc();
    fvm.source() += rPseudoDeltaT*rho.value()*vf.internalField()*mesh().Vsc();

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
dualTimeDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm = scheme_->fvmDdt(rho,vf);
    fvMatrix<Type>& fvm = tfvm.ref();
    const scalarField& rPseudoDeltaT = localRPseudoDeltaT().internalField();

    fvm.diag() += rPseudoDeltaT*rho*mesh().Vsc();
    fvm.source() += rPseudoDeltaT*rho*vf.internalField()*mesh().Vsc();

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
dualTimeDdtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm = scheme_->fvmDdt(rho,vf);
    fvMatrix<Type>& fvm = tfvm.ref();
    const scalarField& rPseudoDeltaT = localRPseudoDeltaT().internalField();

    fvm.diag() += rPseudoDeltaT*alpha*rho*mesh().Vsc();
    fvm.source() += rPseudoDeltaT*alpha*rho*vf.internalField()*
                        mesh().Vsc();

    return tfvm;
}


template<class Type>
tmp<typename dualTimeDdtScheme<Type>::fluxFieldType>
dualTimeDdtScheme<Type>::fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    return scheme_->fvcDdtUfCorr(U,Uf);
}


template<class Type>
tmp<typename dualTimeDdtScheme<Type>::fluxFieldType>
dualTimeDdtScheme<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    return scheme_->fvcDdtPhiCorr(U,phi);
}


template<class Type>
tmp<typename dualTimeDdtScheme<Type>::fluxFieldType>
dualTimeDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    return scheme_->fvcDdtUfCorr(rho,U,Uf);
}


template<class Type>
tmp<typename dualTimeDdtScheme<Type>::fluxFieldType>
dualTimeDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    return scheme_->fvcDdtPhiCorr(rho,U,phi);
}


template<class Type>
tmp<surfaceScalarField> dualTimeDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return scheme_->meshPhi(vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
