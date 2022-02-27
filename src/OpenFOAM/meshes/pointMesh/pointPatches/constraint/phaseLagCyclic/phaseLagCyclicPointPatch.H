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

#ifndef phaseLagCyclicPointPatch_H
#define phaseLagCyclicPointPatch_H

#include "cyclicPointPatch.H"
#include "phaseLagCyclicPolyPatch.H"
#include "pointBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cyclicPointPatch Declaration
\*---------------------------------------------------------------------------*/

class phaseLagCyclicPointPatch
:
    public cyclicPointPatch
{
    // Private data

        //- Local reference cast into the cyclic patch
        const phaseLagCyclicPolyPatch& phaseLagCyclicPolyPatch_;


    // Private Member Functions

        //- No copy construct
        phaseLagCyclicPointPatch(const phaseLagCyclicPointPatch&) = delete;

        //- No copy assignment
        void operator=(const phaseLagCyclicPointPatch&) = delete;


    // Demand driven private data

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName(phaseLagCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        phaseLagCyclicPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    //- Destructor
    virtual ~phaseLagCyclicPointPatch();


    // Member Functions

        // Access

            //- Return the constraint type this pointPatch implements.
            virtual const word& constraintType() const
            {
                return type();
            }

            //- Return the underlying cyclicPolyPatch
            const phaseLagCyclicPolyPatch& phaseLagCyclicPatch() const
            {
                return phaseLagCyclicPolyPatch_;
            }

            //- Return neighbour point patch
            const phaseLagCyclicPointPatch& neighbPatch() const
            {
                label patchi = phaseLagCyclicPolyPatch_.neighbPatchID();
                const pointPatch& pp = this->boundaryMesh()[patchi];
                return refCast<const phaseLagCyclicPointPatch>(pp);
            }

            //- Are the cyclic planes parallel
            bool parallel() const
            {
                return phaseLagCyclicPolyPatch_.parallel();
            }

            //- Return face transformation tensor
            const tensorField& forwardT() const
            {
                return phaseLagCyclicPolyPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            const tensorField& reverseT() const
            {
                return phaseLagCyclicPolyPatch_.reverseT();
            }


        // Access functions for demand driven data

            //- Return the set of pairs of points that require transformation
            //  and/or mapping. First index is on this patch, second on the
            //  neighbour patch.
            virtual const edgeList& transformPairs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
