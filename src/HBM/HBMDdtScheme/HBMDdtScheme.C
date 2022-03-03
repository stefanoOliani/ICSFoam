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

#include "HBMDdtScheme.H"

#include "fvcDiv.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
HBMDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    return tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        IOobject
        (
            "ddt("+dt.name()+')',
            mesh().time().timeName(),
            mesh()
        ),
        mesh(),
        dimensioned<Type>(dt.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
HBMDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        IOobject
        (
            "ddt("+vf.name()+')',
            mesh().time().timeName(),
            mesh()
        ),
        mesh(),
        dimensioned<Type>(vf.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
HBMDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        IOobject
        (
            "ddt("+rho.name()+','+vf.name()+')',
            mesh().time().timeName(),
            mesh()
        ),
        mesh(),
        dimensioned<Type>(rho.dimensions()*vf.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
HBMDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        IOobject
        (
            "ddt("+rho.name()+','+vf.name()+')',
            mesh().time().timeName(),
            mesh()
        ),
        mesh(),
        dimensioned<Type>(rho.dimensions()*vf.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
HBMDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        IOobject
        (
            "ddt("+alpha.name()+','+rho.name()+','+vf.name()+')',
            mesh().time().timeName(),
            mesh()
        ),
        mesh(),
        dimensioned<Type>(rho.dimensions()*vf.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<fvMatrix<Type>>
HBMDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
HBMDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
HBMDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
HBMDdtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            alpha.dimensions()*rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );

    return tfvm;
}


template<class Type>
tmp<typename HBMDdtScheme<Type>::fluxFieldType>
HBMDdtScheme<Type>::fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    tmp<fluxFieldType> tCorr
    (
        new fluxFieldType
        (
            IOobject
            (
                "ddtCorr(" + U.name() + ',' + Uf.name() + ')',
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensioned<typename flux<Type>::type>
            (
                Uf.dimensions()*dimArea/dimTime, Zero
            )
        )
    );

    tCorr.ref().setOriented();

    return tCorr;
}


template<class Type>
tmp<typename HBMDdtScheme<Type>::fluxFieldType>
HBMDdtScheme<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    tmp<fluxFieldType> tCorr
    (
        new fluxFieldType
        (
            IOobject
            (
                "ddtCorr(" + U.name() + ',' + phi.name() + ')',
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensioned<typename flux<Type>::type>
            (
                phi.dimensions()/dimTime, Zero
            )
        )
    );

    tCorr.ref().setOriented();

    return tCorr;
}


template<class Type>
tmp<typename HBMDdtScheme<Type>::fluxFieldType>
HBMDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    tmp<fluxFieldType> tCorr
    (
        new fluxFieldType
        (
            IOobject
            (
                "ddtCorr("
              + rho.name()
              + ',' + U.name() + ',' + Uf.name() + ')',
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensioned<typename flux<Type>::type>
            (
                Uf.dimensions()*dimArea/dimTime, Zero
            )
        )
    );

    tCorr.ref().setOriented();

    return tCorr;
}


template<class Type>
tmp<typename HBMDdtScheme<Type>::fluxFieldType>
HBMDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    tmp<fluxFieldType> tCorr
    (
        new fluxFieldType
        (
            IOobject
            (
                "ddtCorr("
              + rho.name()
              + ',' + U.name() + ',' + phi.name() + ')',
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensioned<typename flux<Type>::type>
            (
                phi.dimensions()/dimTime, Zero
            )
        )
    );

    tCorr.ref().setOriented();

    return tCorr;
}


template<class Type>
tmp<surfaceScalarField> HBMDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return mesh().phi();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
