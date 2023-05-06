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
    Foam::overlapAMIFvPatch

Description
    Adapted from overlapGGI patch implemented in foam-extend by Hrvoje Jasak, 
    Fethi Tekin and Oliver Borm

SourceFiles
    overlapAMIFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef overlapAMIFvPatch_H
#define overlapAMIFvPatch_H

#include "coupledFvPatch.H"
#include "overlapAMILduInterface.H"
#include "overlapAMIPolyPatch.H"
#include "fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class overlapAMIFvPatch Declaration
\*---------------------------------------------------------------------------*/

class overlapAMIFvPatch
:
    public coupledFvPatch,
    public overlapAMILduInterface
{
    // Private data

        const overlapAMIPolyPatch& overlapAMIPolyPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;

        //- Correct patch deltaCoeffs
        virtual void makeDeltaCoeffs(scalarField&) const;

        //- Correct patch non-ortho deltaCoeffs
        virtual void makeNonOrthoDeltaCoeffs(scalarField&) const;

        //- Correct patch non-ortho correction vectors
        virtual void makeNonOrthoCorrVectors(vectorField&) const;

        //- Correct patches after moving points
        virtual void movePoints();


public:

    //- Runtime type information
    TypeName(overlapAMIPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        overlapAMIFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            coupledFvPatch(patch, bm),
            overlapAMILduInterface(),
            overlapAMIPolyPatch_(refCast<const overlapAMIPolyPatch>(patch))
        {}


    // Member functions

        // Access

            //- Return local reference cast into the cyclic patch
            const overlapAMIPolyPatch& overlapAMIPatch() const
            {
                return overlapAMIPolyPatch_;
            }

            //- Return neighbour
            virtual label neighbPatchID() const
            {
                return overlapAMIPolyPatch_.neighbPatchID();
            }

            virtual bool owner() const
            {
                return overlapAMIPolyPatch_.owner();
            }

            //- Return processor number
            virtual const overlapAMIFvPatch& neighbPatch() const
            {
                return refCast<const overlapAMIFvPatch>
                (
                    this->boundaryMesh()[overlapAMIPolyPatch_.neighbPatchID()]
                );
            }

            //- Return a reference to the AMI interpolator
            virtual const AMIPatchToPatchInterpolation& AMI() const
            {
                return overlapAMIPolyPatch_.AMI();
            }

            //- Return true if applying the low weight correction
            virtual bool applyLowWeightCorrection() const
            {
                return overlapAMIPolyPatch_.applyLowWeightCorrection();
            }


            //- Are the cyclic planes parallel
            virtual bool parallel() const
            {
                return overlapAMIPolyPatch_.parallel();
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return overlapAMIPolyPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return overlapAMIPolyPatch_.reverseT();
            }

            const overlapAMIFvPatch& neighbFvPatch() const
            {
                return refCast<const overlapAMIFvPatch>
                (
                    this->boundaryMesh()[overlapAMIPolyPatch_.neighbPatchID()]
                );
            }

            //- Return true if this patch is coupled. This is equivalent
            //- to the coupledPolyPatch::coupled() if parallel running or
            //- both sides present, false otherwise
            virtual bool coupled() const;

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;

            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues = UList<Type>()
            ) const
            {
                return overlapAMIPolyPatch_.interpolate(fld, defaultValues);
            }

            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const tmp<Field<Type>>& tFld,
                const UList<Type>& defaultValues = UList<Type>()
            ) const
            {
                return overlapAMIPolyPatch_.interpolate(tFld, defaultValues);
            }

            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const Field<Type>& fld,
				label cmpt,
                const UList<Type>& defaultValues = UList<Type>()
            ) const
            {
                return overlapAMIPolyPatch_.interpolate(fld, cmpt, defaultValues);
            }

            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const tmp<Field<Type>>& tFld,
				label cmpt,
                const UList<Type>& defaultValues = UList<Type>()
            ) const
            {
                return overlapAMIPolyPatch_.interpolate(tFld, cmpt, defaultValues);
            }

            template<class Type>
			tmp<Field<Type>> untransfInterp
			(
				const Field<Type>& fld,
				const UList<Type>& defaultValues = UList<Type>()
			) const
			{
				return overlapAMIPolyPatch_.untransfInterp(fld, defaultValues);
			}

			template<class Type>
			tmp<Field<Type>> untransfInterp
			(
				const tmp<Field<Type>>& tFld,
				const UList<Type>& defaultValues = UList<Type>()
			) const
			{
				return overlapAMIPolyPatch_.untransfInterp(tFld, defaultValues);
			}

        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //- the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            virtual tmp<labelField> interfaceInternalField
			(
				const labelUList& internalData,
				const labelUList& faceCells
			) const;

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
