/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa

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

Class
    Foam::blockFvMatrix

Description
    Matrix class for block matrices.

SourceFiles
    blockFvMatrix.C

Authors


\*---------------------------------------------------------------------------*/

#ifndef blockFvMatrix_H
#define blockFvMatrix_H

#include "LduMatrix.H"
#include "tmp.H"
#include "fvMatricesFwd.H"
#include "surfaceFields.H"

#include "processorFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class sourceType, class blockType>
class blockFvMatrix;

template<class Type2>
tmp<blockFvMatrix<vector,Type2> > operator*(const tmp<blockFvMatrix<scalar,scalar> >& tA, const Type2& B);


/*---------------------------------------------------------------------------*\
                           Class blockFvMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class sourceType, class blockType>
class blockFvMatrix
:
    public refCount,
    public LduMatrix<sourceType, blockType, blockType>
{
    // Private data

protected:

    // Protected Member Functions

		//- Add patch contribution to internal field
		template<class Type2>
		void addToInternalField
		(
			const labelUList& addr,
			const Field<Type2>& pf,
			Field<Type2>& intf
		) const;

		template<class Type2>
		void addToInternalField
		(
			const labelUList& addr,
			const tmp<Field<Type2>>& tpf,
			Field<Type2>& intf
		) const;

		//- Subtract patch contribution from internal field
		template<class Type2>
		void subtractFromInternalField
		(
			const labelUList& addr,
			const Field<Type2>& pf,
			Field<Type2>& intf
		) const;

		template<class Type2>
		void subtractFromInternalField
		(
			const labelUList& addr,
			const tmp<Field<Type2>>& tpf,
			Field<Type2>& intf
		) const;

		// Matrix completion functionality

			void addBoundaryDiag
			(
				fvMatrix<sourceType> matrix,
				scalarField& diag,
				const direction cmpt
			) const;

			void addBoundarySource
			(
				fvMatrix<sourceType> matrix,
				Field<sourceType>& source,
				const bool couples=true
			) const;

public:

    ClassName("blockFvMatrix");


    // Constructors

        //- Construct given a field to solve for
        blockFvMatrix(const lduMesh&);

        //- Construct as copy
        blockFvMatrix(const blockFvMatrix<sourceType,blockType>&);

        //- Clone
        tmp<blockFvMatrix<sourceType,blockType>> clone() const;

    //- Destructor
    virtual ~blockFvMatrix();

    // Member functions

    	void insertBlock
		(
			const GeometricField<blockType, fvsPatchField, surfaceMesh>& leftField,
			const GeometricField<blockType, fvsPatchField, surfaceMesh>& rightField
		);

    	void insertDissipationBlock
		(
			const GeometricField<blockType, fvsPatchField, surfaceMesh>& dissField
		);

        template<class psiType>
        void Amul
        (
            Field<sourceType>& Apsi,
			const GeometricField<psiType, fvPatchField, volMesh>& psi,
			const fvMesh& mesh
        ) const;

        //Used for block Jacobi smoother
        template<class psiType>
        void AmulNoDiag
        (
            Field<sourceType>& Apsi,
			const GeometricField<psiType, fvPatchField, volMesh>& psi,
			const fvMesh& mesh
        ) const;


        //- Insert fvMatrix into this class (used if this is a
        // diagonal block without inter-component coupling
        // e.g. momentum equation)
        void insertEquation
        (
            fvMatrix<sourceType>& matrix
        );

        template<class psiType>
        void initProcessorInterfaces
		(
			const Field<psiType>& psiif,
			Field<psiType>& result,
			const LduInterfaceFieldPtrsList<psiType>& interfaces
		) const;

        template<class psiType>
        void updateProcessorInterfaces
		(
			Field<sourceType>& result,
			const Pstream::commsTypes commsType,
			const processorFvPatchField<psiType>& procField,
			const fvMesh& mesh,
			const label patchi
		) const;


    //- Operators

        void operator=(const blockFvMatrix<sourceType, blockType>&);
        void operator=(const tmp<blockFvMatrix<sourceType, blockType> >&);
        blockFvMatrix<sourceType, blockType>& operator+=(const blockFvMatrix<sourceType, blockType>&);
        blockFvMatrix<sourceType, blockType>& operator+=(const tmp<blockFvMatrix<sourceType, blockType> >&);
        blockFvMatrix<sourceType, blockType>& operator-=(const blockFvMatrix<sourceType, blockType>&);
        blockFvMatrix<sourceType, blockType>& operator-=(const tmp<blockFvMatrix<sourceType, blockType> >&);

        //- Friends

        template<class Type2>
        friend tmp<blockFvMatrix<vector, Type2> > operator*(const tmp<blockFvMatrix<scalar, scalar> >& tA, const Type2& B);


};


template<class sourceType, class blockType>
tmp<blockFvMatrix<sourceType, blockType> > operator+
(
    const blockFvMatrix<sourceType, blockType>&,
    const blockFvMatrix<sourceType, blockType>&
);

template<class sourceType, class blockType>
tmp<blockFvMatrix<sourceType, blockType> > operator+
(
    const tmp<blockFvMatrix<sourceType, blockType> >&,
    const blockFvMatrix<sourceType, blockType>&
);

template<class sourceType, class blockType>
tmp<blockFvMatrix<sourceType, blockType> > operator+
(
    const blockFvMatrix<sourceType, blockType>&,
    const tmp<blockFvMatrix<sourceType, blockType> >&
);

template<class sourceType, class blockType>
tmp<blockFvMatrix<sourceType, blockType> > operator+
(
    const tmp<blockFvMatrix<sourceType, blockType> >&,
    const tmp<blockFvMatrix<sourceType, blockType> >&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
	#include "blockFvMatrix.C"
	#include "blockFvMatrixUpdateMatrixInterfaces.C"
#endif

// Specialisation for scalars
#include "blockFvScalarScalarMatrix.H"

// Specialisation for vectors
#include "blockFvTensorVectorMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
