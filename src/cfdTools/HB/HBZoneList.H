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
    Foam::HBZoneList

Description
    List container for HB zones

SourceFiles
    HBZoneList.C

\*---------------------------------------------------------------------------*/

#ifndef HBZoneList_H
#define HBZoneList_H

#include "fvMesh.H"
#include "dictionary.H"
#include "fvMatricesFwd.H"
#include "HBZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class HBZoneList;
Ostream& operator<<(Ostream& os, const HBZoneList& models);

/*---------------------------------------------------------------------------*\
                         Class HBZoneList Declaration
\*---------------------------------------------------------------------------*/

class HBZoneList
:
    public PtrList<HBZone>
{
    // Private Member Functions

        //- No copy construct
        HBZoneList(const HBZoneList&) = delete;

        //- No copy assignment
        void operator=(const HBZoneList&) = delete;

        //- Find the period corresponding to the minimum frequence
        scalar findT0() const;

        //- The selected snapshots are common for all the zones
        scalarList selectedSnapshots_;


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;


public:

    //- Constructor
    HBZoneList(const fvMesh& mesh, const dictionary& dict);

    //- Destructor
    ~HBZoneList();


    // Member Functions

        //- Return active status
        bool active(const bool warn = false) const;

        //- Reset the source list
        void reset(const dictionary& dict);

        //- Return const access to the selected snapshots
        const scalarList& selectedSnapshots() const;

        //- Find the set of instants to be used by the HBM
        void setInstants();

        //- Add the HB source term to the specified field
        template<class Type>
        void addingSource
        (
			PtrList<Field<Type>>& source,
			PtrList<PtrList<GeometricField<Type, fvPatchField, volMesh>>>& vars,
			label varNum
        ) const;

        void addingSource
        (
			PtrList<vectorField>& source,
			PtrList<PtrList<volVectorField>>& vars,
			label varNum,
			bool cylCoords
        ) const;

        //- Add the HB source term to the specified equation
        void addingSource
        (
            coupledMatrix& eqSystem,
			PtrList<PtrList<volScalarField>>& scalarVars,
			PtrList<PtrList<volVectorField>>& vectorVars,
			label instantNo
        ) const;

        //- Add the HB source term to the specified equation
        template<class Type>
        void addingSource
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
			bool cylCoords
        ) const;


        //- Add the HB source term to the specified equation
        template<class Type>
        void addingSource
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
        template<class Type>
        void factStep
        (
			PtrList<GeometricField<Type, fvPatchField, volMesh>>& fieldPtr,
			PtrList<scalarField>& deltaTField
        ) const;

        void factStep
        (
			PtrList<PtrList<volScalarField>>& scalarVars,
			PtrList<PtrList<volVectorField>>& vectorVars,
			PtrList<volScalarField>& deltaTField
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
			bool testCylCoord
        ) const;

        // I-O

            //- Read dictionary
            bool read(const dictionary& dict);

            //- Write data to Ostream
            bool writeData(Ostream& os) const;

            //- Ostream operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const HBZoneList& models
            );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HBZoneListTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
