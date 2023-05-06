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

#include "../phaseLagAMIPointPatch/phaseLagAMIPointPatch.H"

#include "pointBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseLagAMIPointPatch, 0);
    addToRunTimeSelectionTable
    (
        facePointPatch,
        phaseLagAMIPointPatch,
        polyPatch
    );
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::phaseLagAMIPointPatch::initGeometry(PstreamBuffers&)
{}


void Foam::phaseLagAMIPointPatch::calcGeometry(PstreamBuffers&)
{}


void Foam::phaseLagAMIPointPatch::initMovePoints
(
    PstreamBuffers&,
    const pointField&
)
{}


void Foam::phaseLagAMIPointPatch::movePoints(PstreamBuffers&, const pointField&)
{}


void Foam::phaseLagAMIPointPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::initUpdateMesh(pBufs);
//    phaseLagAMIPointPatch::initGeometry(pBufs);
}


void Foam::phaseLagAMIPointPatch::updateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::updateMesh(pBufs);
//    phaseLagAMIPointPatch::calcGeometry(pBufs);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseLagAMIPointPatch::phaseLagAMIPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    coupledFacePointPatch(patch, bm),
    phaseLagAMIPolyPatch_(refCast<const phaseLagAMIPolyPatch>(patch))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseLagAMIPointPatch::~phaseLagAMIPointPatch()
{}


// ************************************************************************* //
