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

#ifndef blockFvOperators_H
#define blockFvOperators_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "surfaceInterpolationScheme.H"
#include "blockFvMatrix.H"
#include "tmp.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvm functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvm
{

	tmp<blockFvMatrix<vector,vector>> grad
	(
		const volScalarField& vf,
		const word& name
	);

	tmp<blockFvMatrix<vector,vector>> grad
	(
		const volScalarField& vf
	);

	tmp<blockFvMatrix<scalar,vector>> div
	(
		const volVectorField& vf
	);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Namespace of operators for density-based solvers
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
	sourceType dummyValue  //Needed only to decide the source type
);

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
	sourceType dummyValue  //Needed only to decide the source type
);

template<class blockType>
tmp<blockFvMatrix<blockType, blockType> >
div
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const GeometricField<blockType, fvsPatchField, surfaceMesh>& sf
);

template<class blockType>
tmp<blockFvMatrix<blockType, blockType> >
div
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const tmp<GeometricField<blockType, fvsPatchField, surfaceMesh> >& tsf
);

template<class sourceType, class blockType>
tmp<blockFvMatrix<sourceType, typename outerProduct<blockType,vector>::type> >
grad
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const GeometricField<blockType, fvPatchField, volMesh>& vf,
	sourceType dummyValue  //Needed only to decide the source type
);

template<class sourceType, class blockType>
tmp<blockFvMatrix<sourceType, typename outerProduct<blockType,vector>::type> >
grad
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const tmp<GeometricField<blockType, fvPatchField, volMesh> >& tvf,
	sourceType dummyValue  //Needed only to decide the source type
);


tmp<blockFvMatrix<scalar, vector> >
grad
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const geometricOneField&
);

template<class blockType>
tmp<blockFvMatrix<blockType, typename outerProduct<blockType,vector>::type> >
gradT
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const GeometricField<blockType, fvPatchField, volMesh>& vf
);

template<class blockType>
tmp<blockFvMatrix<blockType, typename outerProduct<blockType,vector>::type> >
gradT
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& w,
    const tmp<GeometricField<blockType, fvPatchField, volMesh> >& tvf
);

template<class blockType>
tmp<blockFvMatrix<blockType, blockType> >
laplacian
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& sf,
    const GeometricField<blockType, fvPatchField, volMesh>& vf,
    bool includePhysicalBoundaries=false
);

template<class blockType>
tmp<blockFvMatrix<blockType, blockType> >
laplacian
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& sf,
    const tmp<GeometricField<blockType, fvPatchField, volMesh> >& tvf,
    bool includePhysicalBoundaries=false
);

template<class blockType>
tmp<blockFvMatrix<blockType, blockType> >
laplacian
(
    const tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> >& tsf,
    const GeometricField<blockType, fvPatchField, volMesh>& vf,
    bool includePhysicalBoundaries=false
);

template<class blockType>
tmp<blockFvMatrix<blockType, blockType> >
laplacian
(
    const tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> >& tsf,
    const tmp<GeometricField<blockType, fvPatchField, volMesh> >& tvf,
    bool includePhysicalBoundaries=false
);

template<class blockType>
tmp<blockFvMatrix<blockType, blockType> >
laplacian
(
    const GeometricField<blockType, fvsPatchField, surfaceMesh>& sf,
    const geometricOneField&,
    bool includePhysicalBoundaries=false
);

template<class blockType>
tmp<blockFvMatrix<blockType, blockType> >
laplacian
(
    const tmp<GeometricField<blockType, fvsPatchField, surfaceMesh> >& tsf,
    const geometricOneField&,
    bool includePhysicalBoundaries=false
);


template<class sourceType, class blockType>
tmp<blockFvMatrix<sourceType, blockType> >
laplacian
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& sf,
    const GeometricField<blockType, fvPatchField, volMesh>& tvf,
	sourceType dummyValue  //Needed only to decide the source type
);


template<class sourceType, class blockType>
tmp<blockFvMatrix<sourceType, blockType> >
laplacian
(
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& sf,
    const tmp<GeometricField<blockType, fvPatchField, volMesh> >& tvf,
	sourceType dummyValue  //Needed only to decide the source type
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvj

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "blockFvOperatorsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
