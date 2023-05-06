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
    Foam::cyclicFvPatch

Description
    Cyclic-plane patch.

SourceFiles
    cyclicFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef phaseLagCyclicFvPatch_H
#define phaseLagCyclicFvPatch_H


#include "phaseLagCyclicPolyPatch.H"
#include "cyclicFvPatch.H"
#include "fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

class phaseLagCyclicFvPatch
:
    public cyclicFvPatch
{
    // Private data

        const phaseLagCyclicPolyPatch& phaseLagCyclicPolyPatch_;


public:

    //- Runtime type information
    TypeName("phaseLagCyclic");


    // Constructors

        //- Construct from polyPatch
        phaseLagCyclicFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            cyclicFvPatch(patch, bm),
            phaseLagCyclicPolyPatch_(refCast<const phaseLagCyclicPolyPatch>(patch))
        {}


    // Member functions

        // Access

            //- Return local reference cast into the cyclic patch
            const phaseLagCyclicPolyPatch& phaseLagCyclicPatch() const
            {
                return phaseLagCyclicPolyPatch_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
