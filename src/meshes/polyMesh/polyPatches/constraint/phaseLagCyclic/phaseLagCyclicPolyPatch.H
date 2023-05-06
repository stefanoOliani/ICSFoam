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
    Foam::phaseLagCyclicPolyPatch

Description
    Cyclic plane patch.

    Note: morph patch face ordering uses geometric matching so with the
    following restrictions:
        -coupled patches should be flat planes.
        -no rotation in patch plane

    Uses coupledPolyPatch::calcFaceTol to calculate
    tolerance per face which might need tweaking.

    Switch on 'phaseLagCyclicPolyPatch' debug flag to write .obj files to show
    the matching.

SourceFiles
    phaseLagCyclicPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef phaseLagCyclicPolyPatch_H
#define phaseLagCyclicPolyPatch_H

#include "cyclicPolyPatch.H"
#include "edgeList.H"
#include "polyBoundaryMesh.H"
#include "diagTensorField.H"
#include "coupleGroupIdentifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class phaseLagCyclicPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class phaseLagCyclicPolyPatch
:
    public cyclicPolyPatch
{
    // Private data

		//- Inter-blade phase angle
		scalar IBPA_;

		//- Name of HB zone
		word HBZoneName_;

		bool cylCoords_;

public:

    //- Runtime type information
    TypeName("phaseLagCyclic");


    // Constructors

        //- Construct from components
        phaseLagCyclicPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const transformType transform = UNKNOWN
        );

        //- Construct from dictionary
        phaseLagCyclicPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        phaseLagCyclicPolyPatch(const phaseLagCyclicPolyPatch&, const polyBoundaryMesh&);

        //- Construct given the original patch and a map
        phaseLagCyclicPolyPatch
        (
            const phaseLagCyclicPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new phaseLagCyclicPolyPatch(*this, bm));
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
                new phaseLagCyclicPolyPatch(*this, bm, index, mapAddressing, newStart)
            );
        }


    //- Destructor
    virtual ~phaseLagCyclicPolyPatch();


    // Member Functions

        const phaseLagCyclicPolyPatch& neighbPatch() const
        {
            const polyPatch& pp = this->boundaryMesh()[this->neighbPatchID()];
            return refCast<const phaseLagCyclicPolyPatch>(pp);
        }

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

        //- Cylindrical coordinates
        const bool& cylCoords() const
        {
            return cylCoords_;
        }

        //- Retrieve the time lag between the two sides
        scalar timeLag() const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
