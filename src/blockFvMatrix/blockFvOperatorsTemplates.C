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

#include "fvCFD.H"
#include "blockFvOperators.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fvj
{

//- Divergence of a cell Jacobian field interpolated to faces with
//  the specified weighting

template<class sourceType, class blockType>
tmp
<
	blockFvMatrix
	<
		sourceType,
		typename innerProduct<blockType,vector>::type
	>
>
div
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const GeometricField<blockType, fvPatchField, volMesh>& vf,
	sourceType dummy
)
{
	typedef typename innerProduct<blockType,vector>::type matrixType;

    const fvMesh& mesh = w.mesh();

    tmp<blockFvMatrix<sourceType, matrixType> > tmx
    (
        new blockFvMatrix<sourceType, matrixType>(mesh)
    );

    blockFvMatrix<sourceType, matrixType>& mx = tmx.ref();

    tmp< surfaceVectorField > lsf = -w*mesh.Sf();
    tmp< surfaceVectorField > usf = mesh.Sf() + lsf();

    Field<matrixType>& upp = mx.upper();
    Field<matrixType>& low = mx.lower();
    Field<matrixType>& diag = mx.diag();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        upp[facei] = (vf[nei]&usf()[facei]);
        low[facei] = (vf[own]&lsf()[facei]);
    }

    mx.negSumDiag();

    mx.interfacesUpper().resize(mesh.boundary().size());
    mx.interfacesLower().resize(mesh.boundary().size());

    forAll(mesh.boundary(), patchi)
    {
        mx.interfacesUpper().set(patchi, vf.boundaryField()[patchi]&usf->boundaryField()[patchi]);
        mx.interfacesLower().set(patchi, vf.boundaryField()[patchi].patchInternalField()&lsf->boundaryField()[patchi]);

        // Don't include physical boundaries because those are dealt with separately
        if (mesh.boundary()[patchi].coupled())
        {
            Field<matrixType>& low =
                mx.interfacesLower()[patchi];

            const labelUList& bOwner =
                mesh.boundary()[patchi].faceCells();

            forAll(bOwner, facei)
            {
                const label own = bOwner[facei];
                diag[own] -= low[facei];
            }
        }
    }
    return tmx;
}

template<class sourceType, class blockType>
tmp
<
	blockFvMatrix
	<
		sourceType,
		typename innerProduct<blockType,vector>::type
	>
>
div
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const tmp<GeometricField<blockType, fvPatchField, volMesh> >& tvf,
	sourceType dummy
)
{
	typedef typename innerProduct<blockType,vector>::type matrixType;

    tmp<blockFvMatrix<sourceType, matrixType> > tmx = div(w, tvf(), dummy);

    tvf.clear();

    return tmx;
}


//- Divergence of a face Jacobian field with specified weightings.
//  Assume surface field already dotted with Sf

template<class blockType>
tmp<blockFvMatrix<blockType, blockType> >
div
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const GeometricField<blockType, fvsPatchField, surfaceMesh>& sf
)
{
    const fvMesh& mesh = w.mesh();

    tmp<blockFvMatrix<blockType, blockType> > tmx
    (
        new blockFvMatrix<blockType, blockType>(mesh)
    );

    blockFvMatrix<blockType, blockType>& mx = tmx.ref();

    Field<blockType>& upp = mx.upper();
    Field<blockType>& low = mx.lower();
    Field<blockType>& diag = mx.diag();

    forAll(upp, facei)
    {
        upp[facei] = sf[facei]*(1-w[facei]);
        low[facei] = -sf[facei]*w[facei];
    }

//    forAll(upp, facei)
//    {
//        upp[facei] = sf[facei]*0.5;
//        low[facei] = -sf[facei]*0.5;
//    }

    mx.negSumDiag();

    mx.interfacesUpper().resize(mesh.boundary().size());
    mx.interfacesLower().resize(mesh.boundary().size());

    forAll(mesh.boundary(), patchi)
    {
        mx.interfacesUpper().set(patchi, sf.boundaryField()[patchi]*(1-w.boundaryField()[patchi]));
        mx.interfacesLower().set(patchi, -sf.boundaryField()[patchi]*w.boundaryField()[patchi]);

//        mx.interfacesUpper().set(patchi, sf.boundaryField()[patchi]*0.5);
//        mx.interfacesLower().set(patchi, -sf.boundaryField()[patchi]*0.5);

        // Don't include physical boundaries because those are dealt with separately
        if (mesh.boundary()[patchi].coupled())
        {
            Field<blockType>& low =
                mx.interfacesLower()[patchi];

            const labelUList& bOwner =
                mesh.boundary()[patchi].faceCells();

            forAll(bOwner, facei)
            {
                const label own = bOwner[facei];
                diag[own] -= low[facei];
            }
        }
    }
    return tmx;
}

template<class blockType>
tmp<blockFvMatrix<blockType, blockType> >
div
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const tmp<GeometricField<blockType, fvsPatchField, surfaceMesh> >& tsf
)
{
    tmp<blockFvMatrix<blockType, blockType> > tmx =
        div(w, tsf());
    tsf.clear();
    return tmx;
}


//- Gradient of a cell Jacobian field interpolated to faces with
//  the specified weighting

template<class sourceType, class blockType>
tmp<blockFvMatrix<sourceType, typename outerProduct<blockType,vector>::type> >
grad
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const GeometricField<blockType, fvPatchField, volMesh>& vf,
	sourceType dummyValue
)
{
    const fvMesh& mesh = w.mesh();

    tmp<blockFvMatrix<sourceType, typename outerProduct<blockType,vector>::type> > tmx
    (
        new blockFvMatrix<sourceType, typename outerProduct<blockType,vector>::type>(mesh)
    );
    blockFvMatrix<sourceType, typename outerProduct<blockType,vector>::type>& mx = tmx.ref();

    tmp< surfaceVectorField > lsf = -w*mesh.Sf();
    tmp< surfaceVectorField > usf = mesh.Sf() + lsf();

    Field<typename outerProduct<blockType,vector>::type>& upp = mx.upper();
    Field<typename outerProduct<blockType,vector>::type>& low = mx.lower();
    Field<typename outerProduct<blockType,vector>::type>& diag = mx.diag();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        upp[facei] = (vf[nei]*usf()[facei]);
        low[facei] = (vf[own]*lsf()[facei]);
    }

    mx.negSumDiag();

    mx.interfacesUpper().resize(mesh.boundary().size());
    mx.interfacesLower().resize(mesh.boundary().size());

    forAll(mesh.boundary(), patchi)
    {
        mx.interfacesUpper().set(patchi, vf.boundaryField()[patchi]*usf->boundaryField()[patchi]);
        mx.interfacesLower().set(patchi, vf.boundaryField()[patchi].patchInternalField()*lsf->boundaryField()[patchi]);

        // Don't include physical boundaries because those are dealt with separately
        if (mesh.boundary()[patchi].coupled())
        {

            Field<typename outerProduct<blockType,vector>::type>& low =
                mx.interfacesLower()[patchi];

            const labelUList& bOwner =
                mesh.boundary()[patchi].faceCells();

            forAll(bOwner, facei)
            {
                const label own = bOwner[facei];
                diag[own] -= low[facei];
            }
        }
    }
    return tmx;
}

template<class sourceType, class blockType>
tmp<blockFvMatrix<sourceType, typename outerProduct<blockType,vector>::type> >
grad
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const tmp<GeometricField<blockType, fvPatchField, volMesh> >& tvf,
	sourceType dummyValue
)
{
    tmp<blockFvMatrix<sourceType, typename outerProduct<blockType,vector>::type> > tmx =
        grad(w, tvf(), dummyValue);

    tvf.clear();

    return tmx;
}


//- Transpose gradient of a cell Jacobian field interpolated to faces with
//- the specified weighting

template<class blockType>
tmp<blockFvMatrix<blockType, typename outerProduct<blockType,vector>::type> >
gradT
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const GeometricField<blockType, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = w.mesh();

    tmp<blockFvMatrix<blockType, typename outerProduct<blockType,vector>::type> > tmx
    (
        new blockFvMatrix<blockType, typename outerProduct<blockType,vector>::type>(mesh)
    );
    blockFvMatrix<blockType, typename outerProduct<blockType,vector>::type>& mx = tmx.ref();

    tmp< surfaceVectorField > lsf = -w*mesh.Sf();
    tmp< surfaceVectorField > usf = mesh.Sf() + lsf();

    Field<typename outerProduct<blockType,vector>::type>& upp = mx.upper();
    Field<typename outerProduct<blockType,vector>::type>& low = mx.lower();
    Field<typename outerProduct<blockType,vector>::type>& diag = mx.diag();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        upp[facei] = (usf()[facei]*vf[nei]);
        low[facei] = (lsf()[facei]*vf[own]);
    }

    mx.negSumDiag();

    mx.interfacesUpper().resize(mesh.boundary().size());
    mx.interfacesLower().resize(mesh.boundary().size());

    forAll(mesh.boundary(), patchi)
    {
        mx.interfacesUpper().set(patchi, usf->boundaryField()[patchi]*vf.boundaryField()[patchi]);
        mx.interfacesLower().set(patchi, lsf->boundaryField()[patchi]*vf.boundaryField()[patchi].patchInternalField());

        // Don't include physical boundaries because those are dealt with separately
        if (mesh.boundary()[patchi].coupled())
        {

            Field<typename outerProduct<blockType,vector>::type>& low =
                mx.interfacesLower()[patchi];

            const labelUList& bOwner =
                mesh.boundary()[patchi].faceCells();

            forAll(bOwner, facei)
            {
                const label own = bOwner[facei];
                diag[own] -= low[facei];
            }
        }
    }
    return tmx;
}

template<class blockType>
tmp<blockFvMatrix<blockType, typename outerProduct<blockType,vector>::type> >
gradT
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const tmp<GeometricField<blockType, fvPatchField, volMesh> >& tvf
)
{
    tmp<blockFvMatrix<blockType, typename outerProduct<blockType,vector>::type> > tmx =
        gradT(w, tvf());
    tvf.clear();
    return tmx;
}


//- Laplacian of a Jacobian volume field

template<class blockType>
tmp<blockFvMatrix<blockType, blockType> >
laplacian
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& sf,
    const GeometricField<blockType, fvPatchField, volMesh>& vf,
    bool includePhysicalBoundaries
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<blockType, blockType> > tmx
    (
        new blockFvMatrix<blockType, blockType>(mesh)
    );
    blockFvMatrix<blockType, blockType>& mx = tmx.ref();

    Field<blockType>& diag = mx.diag();
    Field<blockType>& upp = mx.upper();
    Field<blockType>& low = mx.lower();

    tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > sf2 =
        sf*mesh.magSf()*mesh.deltaCoeffs();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        upp[facei] = vf[nei]*sf2()[facei];
        low[facei] = vf[own]*sf2()[facei];
    }

    mx.negSumDiag();

    mx.interfacesUpper().resize(mesh.boundary().size());
    mx.interfacesLower().resize(mesh.boundary().size());

    forAll(mesh.boundary(), patchi)
    {
        mx.interfacesUpper().set(patchi, vf.boundaryField()[patchi]*sf2().boundaryField()[patchi]);
        mx.interfacesLower().set(patchi, vf.boundaryField()[patchi].patchInternalField()*sf2().boundaryField()[patchi]);

        // Don't include physical boundaries because those are dealt with separately
        if (mesh.boundary()[patchi].coupled() || includePhysicalBoundaries)
        {

            Field<blockType>& low =
                mx.interfacesLower()[patchi];

            const labelUList& bOwner =
                mesh.boundary()[patchi].faceCells();

            forAll(bOwner, facei)
            {
                const label own = bOwner[facei];
                diag[own] -= low[facei];
            }
        }
    }
    return tmx;
}

template<class blockType>
tmp<blockFvMatrix<blockType, blockType> >
laplacian
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& sf,
    const tmp<GeometricField<blockType, fvPatchField, volMesh> >& tvf,
    bool includePhysicalBoundaries
)
{
    tmp<blockFvMatrix<blockType, blockType> > tmx =
        laplacian(sf, tvf(), includePhysicalBoundaries);
    tvf.clear();
    return tmx;
}

template<class blockType>
tmp<blockFvMatrix<blockType, blockType> >
laplacian
(
    const tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> >& tsf,
    const GeometricField<blockType, fvPatchField, volMesh>& vf,
    bool includePhysicalBoundaries
)
{
    tmp<blockFvMatrix<blockType, blockType> > tmx =
        laplacian(tsf(), vf, includePhysicalBoundaries);
    tsf.clear();
    return tmx;
}

template<class blockType>
tmp<blockFvMatrix<blockType, blockType> >
laplacian
(
    const tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> >& tsf,
    const tmp<GeometricField<blockType, fvPatchField, volMesh> >& tvf,
    bool includePhysicalBoundaries
)
{
    tmp<blockFvMatrix<blockType, blockType> > tmx =
        laplacian(tsf(), tvf(), includePhysicalBoundaries);
    tsf.clear();
    tvf.clear();
    return tmx;
}


template<class blockType>
tmp<blockFvMatrix<blockType, blockType> >
laplacian
(
    const GeometricField<blockType, fvsPatchField, surfaceMesh>& sf,
    const geometricOneField&,
    bool includePhysicalBoundaries
)
{
    const fvMesh& mesh = sf.mesh();

    tmp<blockFvMatrix<blockType, blockType> > tmx
    (
        new blockFvMatrix<blockType, blockType>(mesh)
    );
    blockFvMatrix<blockType, blockType>& mx = tmx.ref();

    tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > sf2 =
        sf*mesh.magSf()*mesh.deltaCoeffs();

    Field<blockType>& diag = mx.diag();
    Field<blockType>& upp = mx.upper();
    Field<blockType>& low = mx.lower();

    forAll(upp, facei)
    {
        upp[facei] = sf2()[facei];
        low[facei] = sf2()[facei];  //TODO: don't assign, keep symmetric?
    }

    mx.negSumDiag();

    mx.interfacesUpper().resize(mesh.boundary().size());
    mx.interfacesLower().resize(mesh.boundary().size());

    forAll(mesh.boundary(), patchi)
    {
        mx.interfacesUpper().set(patchi, new Field<blockType>(sf2().boundaryField()[patchi]));
        mx.interfacesLower().set(patchi, new Field<blockType>(sf2().boundaryField()[patchi]));

        // Don't include physical boundaries because those are dealt with explicitly
        if (mesh.boundary()[patchi].coupled() || includePhysicalBoundaries)
        {

            Field<blockType>& low =
                mx.interfacesLower()[patchi];

            const labelUList& bOwner =
                mesh.boundary()[patchi].faceCells();

            forAll(bOwner, facei)
            {
                const label own = bOwner[facei];
                diag[own] -= low[facei];
            }
        }
    }
    return tmx;
}

template<class blockType>
tmp<blockFvMatrix<blockType, blockType> >
laplacian
(
    const tmp<GeometricField<blockType, fvsPatchField, surfaceMesh> >& tsf,
    const geometricOneField&,
    bool includePhysicalBoundaries
)
{
    tmp<blockFvMatrix<blockType, blockType> > tmx =
        laplacian(tsf(), geometricOneField(), includePhysicalBoundaries);
    tsf.clear();
    return tmx;
}


template<class sourceType, class blockType>
tmp<blockFvMatrix<sourceType, blockType> >
laplacian
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& sf,
    const GeometricField<blockType, fvPatchField, volMesh>& vf,
	sourceType dummyValue  //Needed only to decide the source type
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<blockFvMatrix<sourceType, blockType> > tmx
    (
        new blockFvMatrix<sourceType, blockType>(mesh)
    );
    blockFvMatrix<sourceType, blockType>& mx = tmx.ref();

    Field<blockType>& diag = mx.diag();
    Field<blockType>& upp = mx.upper();
    Field<blockType>& low = mx.lower();

    tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> > sf2 =
        sf*mesh.magSf()*mesh.deltaCoeffs();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        upp[facei] = vf[nei]*sf2()[facei];
        low[facei] = vf[own]*sf2()[facei];
    }
    mx.negSumDiag();

    mx.interfacesUpper().resize(mesh.boundary().size());
    mx.interfacesLower().resize(mesh.boundary().size());

    forAll(mesh.boundary(), patchi)
    {
        mx.interfacesUpper().set(patchi, vf.boundaryField()[patchi]*sf2().boundaryField()[patchi]);
        mx.interfacesLower().set(patchi, vf.boundaryField()[patchi].patchInternalField()*sf2().boundaryField()[patchi]);

        // Don't include physical boundaries because those are dealt with separately
        if (mesh.boundary()[patchi].coupled())
        {

            Field<blockType>& low =
                mx.interfacesLower()[patchi];

            const labelUList& bOwner =
                mesh.boundary()[patchi].faceCells();

            forAll(bOwner, facei)
            {
                const label own = bOwner[facei];
                diag[own] -= low[facei];
            }
        }
    }
    return tmx;
}


template<class sourceType, class blockType>
tmp<blockFvMatrix<sourceType, blockType> >
laplacian
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& sf,
    const tmp< GeometricField<blockType, fvPatchField, volMesh> >& tvf,
	sourceType dummyValue  //Needed only to decide the source type
)
{
    tmp<blockFvMatrix<sourceType, blockType> > tmx =
        laplacian(sf, tvf(), dummyValue);

    tvf.clear();

    return tmx;
}


} // End namespace fvj

} // End namespace Foam

// ************************************************************************* //
