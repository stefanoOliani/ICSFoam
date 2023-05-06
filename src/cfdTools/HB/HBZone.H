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

Class
    Foam::HBZone

Description
    HB zone definition based on cell zone and parameters
    obtained from a control dictionary constructed from the given stream.

SourceFiles
    HBZone.C

\*---------------------------------------------------------------------------*/

#ifndef HBZone_H
#define HBZone_H

#include "dictionary.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"
#include "volFieldsFwd.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "mapPolyMesh.H"
#include "Function1.H"
#include "autoPtr.H"
#include "QRMatrix.H"
#include "coupledMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class HBZone Declaration
\*---------------------------------------------------------------------------*/

class HBZone
{
    // Private data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Name of the HB region
        const word name_;

        //- Coefficients dictionary
        dictionary coeffs_;

        //- HB region active flag
        bool active_;

        //- Name of cell zone
        word cellZoneName_;

        //- Cell zone ID
        label cellZoneID_;

        //- Frequency list (rad/s)
        scalarList omegaList_;


        //- Parameters for HB in axisymmetryc cases
        bool cylCoords_;

        vector rotationAxis_;

        point rotationCentre_;

        //- Inverse Fourier transform matrix
        RectangularMatrix<complex> EInv_;

        //- Fourier Matrix (numerically computed from its inverse)
        RectangularMatrix<complex> E_;

        //- HB spectral operator
        scalarSquareMatrix D_;



    // Private Member Functions

        //- No copy construct
        HBZone(const HBZone&) = delete;

        //- No copy assignment
        void operator=(const HBZone&) = delete;

        //- Set the list of frequencies (with their harmonics)
        // for this zone
        void setOmegaList();


public:

    // Declare name of the class and its debug switch
    ClassName("HBZone");


    // Constructors

        //- Construct from fvMesh
        HBZone
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& cellZoneName = word::null
        );

        //- Return clone
        autoPtr<HBZone> clone() const
        {
            NotImplemented;
            return nullptr;
        }


    // Member Functions

        //- Return const access to the HB region name
        inline const word& name() const;

        //- Return const access to the HB active flag
        inline bool active() const;

        //- Return const access to the omega list
        inline const scalarList& omegaList() const;

        inline const vector& rotationAxis() const;

        inline const point& rotationCentre() const;

        //- Return const access to the Fourier matrix
        inline const RectangularMatrix<complex>& E() const;

        //- Return const access to the inverse Fourier matrix
        inline const RectangularMatrix<complex>& EInv() const;

        //- Return const access to the HB operator
        inline const scalarSquareMatrix& D() const;

        //- Calculate the condition number in this zone
        // for the given time instants
        scalar calcConditionN(const scalarList& snapshots);

        //- Update E,EInv and D
        void updateHBOperators(const scalarList& selSnapshots);

        //- Add the HB source term to the specified field
        template<class Type>
        void addSource
        (
			PtrList<Field<Type>>& source,
			PtrList<PtrList<GeometricField<Type, fvPatchField, volMesh>>>& vars,
			label varNum
        ) const;

        void addSource
        (
			PtrList<vectorField>& source,
			PtrList<PtrList<volVectorField>>& vars,
			label varNum,
			bool  tryCylCoords
        ) const;

        //- Add the HB source term to the specified field
        void addSource
        (
			coupledMatrix& eqSystem,
			PtrList<PtrList<volScalarField>>& scalarVars,
			PtrList<PtrList<volVectorField>>& vectorVars,
			label instantNo
        ) const;

        //- Add the HB source term to the specified equation
        template<class Type>
        void addSource
        (
            fvMatrix<Type>& eqn,
			PtrList<GeometricField<Type, fvPatchField, volMesh>>& fieldPtr,
			label instantNo
        ) const;

        //- Add the HB source term to the specified equation
        void addSource
        (
            fvVectorMatrix& eqn,
			PtrList<volVectorField>& fieldPtr,
			label instantNo,
			bool  tryCylCoords
        ) const;

        //- Add the HB source term to the specified equation
        template<class Type>
        void addSource
        (
            const volScalarField& rho,
			fvMatrix<Type>& eqn,
			PtrList<GeometricField<Type, fvPatchField, volMesh>>& fieldPtr,
			label instantNo
        ) const;

        void addBlock
		(
			blockFvMatrix<scalar,scalar>& blockMatrix,
			label rowIndex,
			label colIndex
		) const;

        void addBlock
		(
			blockFvMatrix<vector,tensor>& blockMatrix,
			label rowIndex,
			label colIndex
		) const;

        //- Factorization step to stabilize explicit HB method
        void factorizationStep
        (
			PtrList<PtrList<volScalarField>>& scalarVars,
			PtrList<PtrList<volVectorField>>& vectorVars,
			PtrList<volScalarField>& deltaTField
        ) const;

        //- Factorization step to stabilize explicit HB method
        template<class Type>
        void factorizationStep
        (
			PtrList<GeometricField<Type, fvPatchField, volMesh>>& fieldPtr,
			PtrList<scalarField>& deltaTField
        ) const;

        //- Reconstruct the field at desired time from the snapshot
		template<class Type>
		void reconstruct
		(
			GeometricField<Type, fvPatchField, volMesh>& reconstFld,
			PtrList<GeometricField<Type, fvPatchField, volMesh>>& fieldPtr,
			scalar actTime
		) const;

		void reconstruct
		(
			volVectorField& reconstFld,
			PtrList<volVectorField>& fieldPtr,
			scalar actTime,
			bool tryCylCoords
		) const;

    // I-O

        //- Write
        void writeData(Ostream& os) const;

        //- Read MRF dictionary
        bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HBZoneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "HBZoneI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
