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
    Foam::phaseLagAMIFvPatch

Description
    Phase Lag Cyclic patch for Arbitrary Mesh Interface (AMI)

SourceFiles
    phaseLagAMIFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef phaseLagAMIFvPatch_H
#define phaseLagAMIFvPatch_H

#include "coupledFvPatch.H"
#include "phaseLagAMILduInterface.H"
#include "phaseLagAMIPolyPatch.H"
#include "fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class phaseLagAMIFvPatch Declaration
\*---------------------------------------------------------------------------*/

class phaseLagAMIFvPatch
:
    public coupledFvPatch,
    public phaseLagAMILduInterface
{
    // Private data

        const phaseLagAMIPolyPatch& phaseLagAMIPolyPatch_;


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
    TypeName(phaseLagAMIPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        phaseLagAMIFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            coupledFvPatch(patch, bm),
            phaseLagAMILduInterface(),
            phaseLagAMIPolyPatch_(refCast<const phaseLagAMIPolyPatch>(patch))
        {}


    // Member functions

        // Access

            //- Return local reference cast into the cyclic patch
            const phaseLagAMIPolyPatch& phaseLagAMIPatch() const
            {
                return phaseLagAMIPolyPatch_;
            }

            //- Return neighbour
            virtual label neighbPatchID() const
            {
                return phaseLagAMIPolyPatch_.neighbPatchID();
            }

            virtual bool owner() const
            {
                return phaseLagAMIPolyPatch_.owner();
            }

            //- Return processor number
            virtual const phaseLagAMIFvPatch& neighbPatch() const
            {
                return refCast<const phaseLagAMIFvPatch>
                (
                    this->boundaryMesh()[phaseLagAMIPolyPatch_.neighbPatchID()]
                );
            }

            //- Return a reference to the AMI interpolator
            virtual const AMIPatchToPatchInterpolation& AMI() const
            {
                return phaseLagAMIPolyPatch_.AMI();
            }

            //- Return true if applying the low weight correction
            virtual bool applyLowWeightCorrection() const
            {
                return phaseLagAMIPolyPatch_.applyLowWeightCorrection();
            }


            //- Are the cyclic planes parallel
            virtual bool parallel() const
            {
                return phaseLagAMIPolyPatch_.parallel();
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return phaseLagAMIPolyPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return phaseLagAMIPolyPatch_.reverseT();
            }

            const phaseLagAMIFvPatch& neighbFvPatch() const
            {
                return refCast<const phaseLagAMIFvPatch>
                (
                    this->boundaryMesh()[phaseLagAMIPolyPatch_.neighbPatchID()]
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
				const word& fieldName,
                const UList<Type>& defaultValues = UList<Type>()
            ) const
            {
                return phaseLagAMIPolyPatch_.interpolate(fld, fieldName, defaultValues);
            }

            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const tmp<Field<Type>>& tFld,
				const word& fieldName,
                const UList<Type>& defaultValues = UList<Type>()
            ) const
            {
                return phaseLagAMIPolyPatch_.interpolate(tFld, fieldName, defaultValues);
            }

            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues = UList<Type>()
            ) const
            {
                return phaseLagAMIPolyPatch_.interpolate(fld, defaultValues);
            }

            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const tmp<Field<Type>>& tFld,
                const UList<Type>& defaultValues = UList<Type>()
            ) const
            {
                return phaseLagAMIPolyPatch_.interpolate(tFld, defaultValues);
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
				const labelUList& faceCell
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
