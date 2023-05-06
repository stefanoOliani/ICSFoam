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

\*---------------------------------------------------------------------------*/

#include "../overlapAMIPointPatch/overlapAMIPointPatch.H"

#include "pointBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(overlapAMIPointPatch, 0);
    addToRunTimeSelectionTable
    (
        facePointPatch,
        overlapAMIPointPatch,
        polyPatch
    );
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::overlapAMIPointPatch::initGeometry(PstreamBuffers&)
{}


void Foam::overlapAMIPointPatch::calcGeometry(PstreamBuffers&)
{}


void Foam::overlapAMIPointPatch::initMovePoints
(
    PstreamBuffers&,
    const pointField&
)
{}


void Foam::overlapAMIPointPatch::movePoints(PstreamBuffers&, const pointField&)
{}


void Foam::overlapAMIPointPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::initUpdateMesh(pBufs);
//    overlapAMIPointPatch::initGeometry(pBufs);
}


void Foam::overlapAMIPointPatch::updateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::updateMesh(pBufs);
//    overlapAMIPointPatch::calcGeometry(pBufs);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::overlapAMIPointPatch::overlapAMIPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    coupledFacePointPatch(patch, bm),
    overlapAMIPolyPatch_(refCast<const overlapAMIPolyPatch>(patch))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::overlapAMIPointPatch::~overlapAMIPointPatch()
{}


// ************************************************************************* //
