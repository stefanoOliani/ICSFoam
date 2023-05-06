/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

    Hrvoje Jasak, Wikki Ltd.  All rights reserved
    Fethi Tekin, All rights reserved.
    Oliver Borm, All rights reserved.

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
    Foam::overlapAMIPolyPatch

Description
    Adapted from overlapGGI patch implemented in foam-extend by Hrvoje Jasak, 
    Fethi Tekin and Oliver Borm

SourceFiles
    overlapAMIPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef overlapAMIPolyPatch_H
#define overlapAMIPolyPatch_H

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
                    Class overlapAMIPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class overlapAMIPolyPatch
:
    public coupledPolyPatch
{
    // Private data

	 	 //- Rotation axis
		const vector rotationAxis_;

		// Number of copies in order to complete 360 degrees
		const label nCopies_;

		// Pointer to expanded patch
		mutable autoPtr<primitivePatch> expandedPatchPtr_;

		mutable faceList expandedFaces_;
		mutable pointField expandedPoints_;

    // Private Member Functions

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
    TypeName("overlapAMI");


    // Constructors

        //- Construct from (base coupled patch) components
        overlapAMIPolyPatch
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
        overlapAMIPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const word& defaultAMIMethod = faceAreaWeightAMI::typeName
        );

        //- Construct as copy, resetting the boundary mesh
        overlapAMIPolyPatch(const overlapAMIPolyPatch&, const polyBoundaryMesh&);

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        overlapAMIPolyPatch
        (
            const overlapAMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& nbrPatchName
        );

        //- Construct given the original patch and a map
        overlapAMIPolyPatch
        (
            const overlapAMIPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );


        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>::NewFrom<overlapAMIPolyPatch>(*this, bm);
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
                new overlapAMIPolyPatch
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
                new overlapAMIPolyPatch
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
    virtual ~overlapAMIPolyPatch() = default;


    // Member Functions

        // Access

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

            inline const vector& rotationAxis() const;

			//- Return wedge angle
			inline scalar angle() const;

			//- Return number of slave copies
			inline const label& nCopies() const;

            //- Neighbour patch name
            inline const word& neighbPatchName() const;

            //- Neighbour patch ID
            virtual label neighbPatchID() const;

            //- Does this side own the patch?
            virtual bool owner() const;

            //- Return a reference to the neighbour patch
            virtual const overlapAMIPolyPatch& neighbPatch() const;

            //- Return a reference to the AMI interpolator
            const AMIPatchToPatchInterpolation& AMI() const;

            //- Return a reference to the expanded patch
			const primitivePatch& expandedPatch() const;

            //- Helper function to return the weights
            inline const scalarListList& weights() const;

            //- Helper function to return the weights sum
            inline const scalarField& weightsSum() const;

            //- Return true if applying the low weight correction
            bool applyLowWeightCorrection() const;


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

            // Interpolations

                //- Expand face field for 360 degrees coverage
				template<class Type>
				tmp<Field<Type> > expandData(const Field<Type>& spf) const;

				//- Expand face field for 360 degrees coverage
				template<class Type>
				tmp<Field<Type> > expandData(const Field<Type>& spf, label cmpt) const;

				//- Expand face field for 360 degrees coverage
				template<class Type>
				tmp<Field<Type> > untransfExpandData(const Field<Type>& spf) const;

				template<class Type>
				UList<Type> expandData(const UList<Type>& defValues) const;

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


                //- Interpolate field
                template<class Type>
                tmp<Field<Type>> interpolate
                (
                    const Field<Type>& fld,
					label cmpt,
                    const UList<Type>& defaultValues = UList<Type>()
                ) const;

                //- Interpolate tmp field
                template<class Type>
                tmp<Field<Type>> interpolate
                (
                    const tmp<Field<Type>>& tFld,
					label cmpt,
                    const UList<Type>& defaultValues = UList<Type>()
                ) const;

                //- Low-level interpolate List
                template<class Type, class CombineOp>
                void interpolate
                (
                    const UList<Type>& fld,
					label cmpt,
                    const CombineOp& cop,
                    List<Type>& result,
                    const UList<Type>& defaultValues = UList<Type>()
                ) const;


                //- Interpolate field
			   template<class Type>
			   tmp<Field<Type>> untransfInterp
			   (
				   const Field<Type>& fld,
				   const UList<Type>& defaultValues = UList<Type>()
			   ) const;

			   //- Interpolate tmp field
			   template<class Type>
			   tmp<Field<Type>> untransfInterp
			   (
				   const tmp<Field<Type>>& tFld,
				   const UList<Type>& defaultValues = UList<Type>()
			   ) const;

			   //- Low-level interpolate List
			   template<class Type, class CombineOp>
			   void untransfInterp
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

        // Compute rotation tensor using Rodrigues' Rotation Formula
        // Requires: rotation axis + rotation angle
        tensor RodriguesRotation
        (
            const vector& rotationAxis,
            const scalar& rotationAngle,
            const bool inDegrees = true
        ) const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};

template<>
tmp<scalarField > overlapAMIPolyPatch::expandData(const scalarField& spf, label cmpt) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "overlapAMIPolyPatchI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "overlapAMIPolyPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
