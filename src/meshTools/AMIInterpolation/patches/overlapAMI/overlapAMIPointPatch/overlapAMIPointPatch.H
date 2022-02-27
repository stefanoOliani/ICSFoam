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

#ifndef overlapAMIPointPatch_H
#define overlapAMIPointPatch_H

#include "overlapAMIPolyPatch.H"
#include "coupledFacePointPatch.H"
#include "pointBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class overlapAMIPointPatch Declaration
\*---------------------------------------------------------------------------*/

class overlapAMIPointPatch
:
    public coupledFacePointPatch
{
    // Private data

        //- Local reference cast into the cyclic AMI patch
        const overlapAMIPolyPatch& overlapAMIPolyPatch_;


    // Private Member Functions

        //- No copy construct
        overlapAMIPointPatch(const overlapAMIPointPatch&) = delete;

        //- No copy assignment
        void operator=(const overlapAMIPointPatch&) = delete;


protected:

    // Protected Member Functions

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
    TypeName(overlapAMIPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        overlapAMIPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    //- Destructor
    virtual ~overlapAMIPointPatch();


    // Member Functions

        //- Return the constraint type this pointPatch implements.
        virtual const word& constraintType() const
        {
            return type();
        }

        //- Return the underlying overlapAMIPolyPatch
        const overlapAMIPolyPatch& overlapAMIPatch() const
        {
            return overlapAMIPolyPatch_;
        }

        //- Return neighbour point patch
        const overlapAMIPointPatch& neighbPatch() const
        {
            label patchi = overlapAMIPolyPatch_.neighbPatchID();
            const pointPatch& pp = this->boundaryMesh()[patchi];
            return refCast<const overlapAMIPointPatch>(pp);
        }

        //- Are the cyclic planes parallel
        bool parallel() const
        {
            return overlapAMIPolyPatch_.parallel();
        }

        //- Return face transformation tensor
        const tensorField& forwardT() const
        {
            return overlapAMIPolyPatch_.forwardT();
        }

        //- Return neighbour-cell transformation tensor
        const tensorField& reverseT() const
        {
            return overlapAMIPolyPatch_.reverseT();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
