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

\*---------------------------------------------------------------------------*/

#include "phaseLagCyclicPointPatch.H"
#include "pointBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "edgeList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseLagCyclicPointPatch, 0);
    addToRunTimeSelectionTable
    (
        facePointPatch,
		phaseLagCyclicPointPatch,
        polyPatch
    );
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::phaseLagCyclicPointPatch::initGeometry(PstreamBuffers&)
{}


void Foam::phaseLagCyclicPointPatch::calcGeometry(PstreamBuffers&)
{}


void Foam::phaseLagCyclicPointPatch::initMovePoints(PstreamBuffers&, const pointField&)
{}


void Foam::phaseLagCyclicPointPatch::movePoints(PstreamBuffers&, const pointField&)
{}


void Foam::phaseLagCyclicPointPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::initUpdateMesh(pBufs);
    phaseLagCyclicPointPatch::initGeometry(pBufs);
}


void Foam::phaseLagCyclicPointPatch::updateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::updateMesh(pBufs);
    phaseLagCyclicPointPatch::calcGeometry(pBufs);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseLagCyclicPointPatch::phaseLagCyclicPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    cyclicPointPatch(patch, bm),
	phaseLagCyclicPolyPatch_(refCast<const phaseLagCyclicPolyPatch>(patch))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseLagCyclicPointPatch::~phaseLagCyclicPointPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::edgeList& Foam::phaseLagCyclicPointPatch::transformPairs() const
{
    return phaseLagCyclicPolyPatch_.coupledPoints();
}


// ************************************************************************* //
