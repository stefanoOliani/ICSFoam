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

#ifndef phaseLagAMIPointPatch_H
#define phaseLagAMIPointPatch_H

#include "phaseLagAMIPolyPatch.H"
#include "coupledFacePointPatch.H"
#include "pointBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class phaseLagAMIPointPatch Declaration
\*---------------------------------------------------------------------------*/

class phaseLagAMIPointPatch
:
    public coupledFacePointPatch
{
    // Private data

        //- Local reference cast into the cyclic AMI patch
        const phaseLagAMIPolyPatch& phaseLagAMIPolyPatch_;


    // Private Member Functions

        //- No copy construct
        phaseLagAMIPointPatch(const phaseLagAMIPointPatch&) = delete;

        //- No copy assignment
        void operator=(const phaseLagAMIPointPatch&) = delete;


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
    TypeName(phaseLagAMIPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        phaseLagAMIPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    //- Destructor
    virtual ~phaseLagAMIPointPatch();


    // Member Functions

        //- Return the constraint type this pointPatch implements.
        virtual const word& constraintType() const
        {
            return type();
        }

        //- Return the underlying phaseLagAMIPolyPatch
        const phaseLagAMIPolyPatch& phaseLagAMIPatch() const
        {
            return phaseLagAMIPolyPatch_;
        }

        //- Return neighbour point patch
        const phaseLagAMIPointPatch& neighbPatch() const
        {
            label patchi = phaseLagAMIPolyPatch_.neighbPatchID();
            const pointPatch& pp = this->boundaryMesh()[patchi];
            return refCast<const phaseLagAMIPointPatch>(pp);
        }

        //- Are the cyclic planes parallel
        bool parallel() const
        {
            return phaseLagAMIPolyPatch_.parallel();
        }

        //- Return face transformation tensor
        const tensorField& forwardT() const
        {
            return phaseLagAMIPolyPatch_.forwardT();
        }

        //- Return neighbour-cell transformation tensor
        const tensorField& reverseT() const
        {
            return phaseLagAMIPolyPatch_.reverseT();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
