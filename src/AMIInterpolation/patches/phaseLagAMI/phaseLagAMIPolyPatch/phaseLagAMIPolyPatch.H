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
    Foam::phaseLagAMIPolyPatch

Description
	No topology change, only rotational transform

SourceFiles
    phaseLagAMIPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef phaseLagAMIPolyPatch_H
#define phaseLagAMIPolyPatch_H

#include "coupledPolyPatch.H"
#include "AMIPatchToPatchInterpolation.H"
#include "polyBoundaryMesh.H"
#include "coupleGroupIdentifier.H"
#include "faceAreaWeightAMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class OBJstream;

/*---------------------------------------------------------------------------*\
                    Class phaseLagAMIPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class phaseLagAMIPolyPatch
:
    public coupledPolyPatch
{
    // Private data

        //- Periodic patch name
        mutable word periodicPatchName_;

        //- Periodic patch ID
        mutable label periodicPatchID_;

        word HBZoneName_;

        //- Sub time level number
        mutable label subTimeLevel_;

        //- Inter-blade phase angle
        scalar IBPA_;

        bool cylCoords_;

        //- Current number of transformations
        mutable label nTransformsFwd_;
        mutable label nTransformsBwd_;
        mutable label nTransforms_;

		// Pointer to expanded patch
		mutable autoPtr<primitivePatch> expandedPatchPtr_;

		mutable faceList expandedFaces_;
		mutable pointField expandedPoints_;

        //- Number of sectors in a rotationally periodic geometry (optional)
        const label nSectors_;

        //- Maximum number of attempts to match the AMI geometry
        const label maxIter_;

    // Private Member Functions

        //- Synchronise the periodic transformations
        void syncTransforms() const;

        //- Debug: write obj files of patch (collected on master)
        void writeOBJ(const primitivePatch& p, OBJstream& str) const;

protected:

    // Protected data

        //- Name of other half
        mutable word nbrPatchName_;

        //- Optional patchGroup to find neighbPatch
        const coupleGroupIdentifier coupleGroup_;

        //- Index of other half
        mutable label nbrPatchID_;

        //- AMI interpolation class
        mutable autoPtr<AMIPatchToPatchInterpolation> AMIPtr_;

        //- Dictionary used during projection surface construction
        const dictionary surfDict_;

        //- Projection surface
        mutable autoPtr<searchableSurface> surfPtr_;



    // Protected Member Functions

        //- Create and return pointer to the projection surface
        const autoPtr<searchableSurface>& surfPtr() const;

        //-  Reset the AMI interpolator, use current patch points
        virtual void resetAMI() const;

        //- Recalculate the transformation tensors
        virtual void calcTransforms();

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers& pBufs, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers& pBufs, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);

        //- Clear geometry
        virtual void clearGeom();


public:

    //- Runtime type information
    TypeName("phaseLagAMI");


    // Constructors

        //- Construct from (base coupled patch) components
        phaseLagAMIPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const transformType transform = UNKNOWN,
            const word& defaultAMIMethod = faceAreaWeightAMI::typeName
        );

        //- Construct from dictionary
        phaseLagAMIPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const word& defaultAMIMethod = faceAreaWeightAMI::typeName
        );

        //- Construct as copy, resetting the boundary mesh
        phaseLagAMIPolyPatch(const phaseLagAMIPolyPatch&, const polyBoundaryMesh&);

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        phaseLagAMIPolyPatch
        (
            const phaseLagAMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& nbrPatchName
        );

        //- Construct given the original patch and a map
        phaseLagAMIPolyPatch
        (
            const phaseLagAMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );


        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>::NewFrom<phaseLagAMIPolyPatch>(*this, bm);
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new phaseLagAMIPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    nbrPatchName_
                )
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new phaseLagAMIPolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~phaseLagAMIPolyPatch() = default;


    // Member Functions

        // Access

			//- Periodic patch ID
			virtual label periodicPatchID() const;

	        //- Axis of rotation for rotational phaseLagCyclics
	        const scalar& IBPA() const
	        {
	            return IBPA_;
	        }

	        //- Harmonic balance zone for phaseLagCyclics
	        const word& HBZoneName() const
	        {
	            return HBZoneName_;
	        }

            //- Tolerance used e.g. for area calculations/limits
            static const scalar tolerance_;

            //- Flag to indicate whether the AMI can be reset
            inline bool canResetAMI() const;

            //- Is patch 'coupled'. Note that on AMI the geometry is not
            //- coupled but the fields are!
            virtual bool coupled() const
            {
                return false;
            }

            //- Neighbour patch name
            inline const word& neighbPatchName() const;

            //- Neighbour patch name
			inline label& nTransformsFwd() const;

			inline label& nTransformsBwd() const;

			inline label& nTransforms() const;

            //- Neighbour patch ID
            virtual label neighbPatchID() const;

            //- Does this side own the patch?
            virtual bool owner() const;

            //- Return a reference to the neighbour patch
            virtual const phaseLagAMIPolyPatch& neighbPatch() const;

            //- Return a reference to the AMI interpolator
            const AMIPatchToPatchInterpolation& AMI() const;

            //- Helper function to return the weights
            inline const scalarListList& weights() const;

            //- Helper function to return the weights sum
            inline const scalarField& weightsSum() const;

            //- Return true if applying the low weight correction
            bool applyLowWeightCorrection() const;

            //- Retrieve the time lag between the two sides
            scalar timeLag(bool fwd) const;


            // Transformations

                //- Transform patch-based positions from nbr side to this side
                virtual void transformPosition(pointField&) const;

                //- Transform a patch-based position from nbr side to this side
                virtual void transformPosition
                (
                    point& l,
                    const label facei
                ) const;

                //- Transform a patch-based position from this side to nbr side
                virtual void reverseTransformPosition
                (
                    point& l,
                    const label facei
                ) const;

                //- Transform a patch-based direction from this side to nbr side
                virtual void reverseTransformDirection
                (
                    vector& d,
                    const label facei
                ) const;

                //- Transform patch-based positions from nbr side to this side
				void revTransformPositionPatch
				(
					pointField& l,
					const coupledPolyPatch& perioPatch
				) const;


            // Interpolations

                //- Expand face field for 360 degrees coverage
				template<class Type>
				tmp<Field<Type> > expandData(const Field<Type>& spf) const;

				template<class Type>
				tmp<Field<Type> > expandData(const Field<Type>& spf, const word& fieldName) const;

				template<class Type>
				UList<Type> expandData(const UList<Type>& defValues) const;

                //- Interpolate field
                template<class Type>
                tmp<Field<Type>> interpolate
                (
                    const Field<Type>& fld,
					const word& fieldName,
                    const UList<Type>& defaultValues = UList<Type>()
                ) const;

                //- Interpolate tmp field
                template<class Type>
                tmp<Field<Type>> interpolate
                (
                    const tmp<Field<Type>>& tFld,
					const word& fieldName,
                    const UList<Type>& defaultValues = UList<Type>()
                ) const;

                //- Low-level interpolate List
                template<class Type, class CombineOp>
                void interpolate
                (
                    const UList<Type>& fld,
					const word& fieldName,
                    const CombineOp& cop,
                    List<Type>& result,
                    const UList<Type>& defaultValues = UList<Type>()
                ) const;

                //- Interpolate field
                template<class Type>
                tmp<Field<Type>> interpolate
                (
                    const Field<Type>& fld,
                    const UList<Type>& defaultValues = UList<Type>()
                ) const;

                //- Interpolate tmp field
                template<class Type>
                tmp<Field<Type>> interpolate
                (
                    const tmp<Field<Type>>& tFld,
                    const UList<Type>& defaultValues = UList<Type>()
                ) const;

                //- Low-level interpolate List
                template<class Type, class CombineOp>
                void interpolate
                (
                    const UList<Type>& fld,
                    const CombineOp& cop,
                    List<Type>& result,
                    const UList<Type>& defaultValues = UList<Type>()
                ) const;


        //- Calculate the patch geometry
        virtual void calcGeometry
        (
            const primitivePatch& referPatch,
            const pointField& thisCtrs,
            const vectorField& thisAreas,
            const pointField& thisCc,
            const pointField& nbrCtrs,
            const vectorField& nbrAreas,
            const pointField& nbrCc
        );

        //- Initialize ordering for primitivePatch. Does not
        //- refer to *this (except for name() and type() etc.)
        virtual void initOrder
        (
            PstreamBuffers&,
            const primitivePatch&
        ) const;

        //- Return new ordering for primitivePatch.
        //  Ordering is -faceMap: for every face
        //  index of the new face -rotation:for every new face the clockwise
        //  shift of the original face. Return false if nothing changes
        //  (faceMap is identity, rotation is 0), true otherwise.
        virtual bool order
        (
            PstreamBuffers&,
            const primitivePatch&,
            labelList& faceMap,
            labelList& rotation
        ) const;

        //- Return face index on neighbour patch which shares point p
        //- following trajectory vector n
        label pointFace
        (
            const label facei,
            const vector& n,
            point& p
        ) const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


template<>
tmp<scalarField > phaseLagAMIPolyPatch::expandData(const scalarField& spf, const word& fieldName) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseLagAMIPolyPatchI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "phaseLagAMIPolyPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
