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

#include "overlapAMIPolyPatch.H"
#include "faceAreaWeightAMI.H"
#include "SubField.H"
#include "Time.H"
#include "unitConversion.H"
#include "OFstream.H"
#include "meshTools.H"
#include "addToRunTimeSelectionTable.H"

// For debugging
#include "OBJstream.H"
#include "PatchTools.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(overlapAMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, overlapAMIPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, overlapAMIPolyPatch, dictionary);
}

const Foam::scalar Foam::overlapAMIPolyPatch::tolerance_ = 1e-10;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::overlapAMIPolyPatch::writeOBJ
(
    const primitivePatch& p,
    OBJstream& str
) const
{
    // Collect faces and points
    pointField allPoints;
    faceList allFaces;
    labelList pointMergeMap;
    PatchTools::gatherAndMerge
    (
        -1.0,           // do not merge points
        p,
        allPoints,
        allFaces,
        pointMergeMap
    );

    if (Pstream::master())
    {
        // Write base geometry
        str.write(allFaces, allPoints);
    }
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

const Foam::autoPtr<Foam::searchableSurface>&
Foam::overlapAMIPolyPatch::surfPtr() const
{
    const word surfType(surfDict_.getOrDefault<word>("type", "none"));

    if (!surfPtr_ && owner() && surfType != "none")
    {
        word surfName(surfDict_.getOrDefault("name", name()));

        const polyMesh& mesh = boundaryMesh().mesh();

        surfPtr_ =
            searchableSurface::New
            (
                surfType,
                IOobject
                (
                    surfName,
                    mesh.time().constant(),
                    "triSurface",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                surfDict_
            );
    }

    return surfPtr_;
}


void Foam::overlapAMIPolyPatch::resetAMI() const
{
    if (owner())
    {
    	const scalar myAngleThis = 360.0/scalar(nCopies());
    	const scalar myAngleNbr = 360.0/scalar(neighbPatch().nCopies());

        const pointField& thisPoints0 = localPoints();
        const pointField& nbrPoints0 = neighbPatch().localPoints();

		expandedPoints_.resize(nCopies()*thisPoints0.size());
		expandedFaces_.resize(nCopies()*localFaces().size());

		neighbPatch().expandedPoints_.resize(neighbPatch().nCopies()*nbrPoints0.size());
		neighbPatch().expandedFaces_.resize(neighbPatch().nCopies()*neighbPatch().localFaces().size());

        autoPtr<OBJstream> ownStr;
        autoPtr<OBJstream> neiStr;
        if (debug)
        {
            const Time& runTime = boundaryMesh().mesh().time();

            fileName dir(runTime.globalPath());
            fileName postfix("_" + runTime.timeName()+"_expanded.obj");

            ownStr.reset(new OBJstream(dir/name() + postfix));
            neiStr.reset(new OBJstream(dir/neighbPatch().name() + postfix));

            InfoInFunction
                << "patch:" << name()
                << " writing accumulated AMI to " << ownStr().name()
                << " and " << neiStr().name() << endl;
        }

        //Expand geometry

        //TGT step

        label nPointsGeomThis = 0;
        label nFacesGeomThis = 0;

		for (label i = 0; i < nCopies(); ++ i)
		{
			faceList geomLocalFacesSrc(localFaces());

			// Calculate transform
			const tensor curRotation = RodriguesRotation(rotationAxis_,  i*myAngleThis);

			forAll (thisPoints0, pointI)
			{
				expandedPoints_[nPointsGeomThis] = Foam::transform(curRotation, thisPoints0[pointI]);
				nPointsGeomThis++;
			}

			const label copyOffsetGeomSrc = i*thisPoints0.size();

			forAll (geomLocalFacesSrc, faceI)
			{
				const face& curGeomFace = geomLocalFacesSrc[faceI];

				face& curExpandedFace = expandedFaces_[nFacesGeomThis];

				// Copy face with offsets
				curExpandedFace.setSize(curGeomFace.size());

				forAll (curGeomFace, fpI)
				{
					curExpandedFace[fpI] = curGeomFace[fpI] + copyOffsetGeomSrc;
				}

				nFacesGeomThis++;
			}
		}

		//SRC step

		label nPointsGeomNbr = 0;
		label nFacesGeomNbr = 0;

		for (label i = 0; i < neighbPatch().nCopies(); ++ i)
		{
			faceList geomLocalFacesTgt(neighbPatch().localFaces());

			// Calculate transform
			const tensor curRotation = RodriguesRotation(rotationAxis_,  i*myAngleNbr);

			forAll (nbrPoints0, pointI)
			{
				neighbPatch().expandedPoints_[nPointsGeomNbr] = Foam::transform(curRotation, nbrPoints0[pointI]);
				nPointsGeomNbr++;
			}

			const label copyOffsetGeomNbr = i*nbrPoints0.size();

			forAll (geomLocalFacesTgt, faceI)
			{
				const face& curGeomFace = geomLocalFacesTgt[faceI];

				face& curExpandedFace = neighbPatch().expandedFaces_[nFacesGeomNbr];

				// Copy face with offsets
				curExpandedFace.setSize(curGeomFace.size());

				forAll (curGeomFace, fpI)
				{
					curExpandedFace[fpI] = curGeomFace[fpI] + copyOffsetGeomNbr;
				}

				nFacesGeomNbr++;
			}
		}

        // Close debug streams
        if (ownStr.valid())
        {
            ownStr.clear();
        }
        if (neiStr.valid())
        {
            neiStr.clear();
        }

		SubList<face> expandedFacesSrc(expandedFaces_);
		expandedPatchPtr_.set(new primitivePatch(expandedFacesSrc, expandedPoints_));

		SubList<face> expandedFacesTrg(neighbPatch().expandedFaces_);
		neighbPatch().expandedPatchPtr_.set(new primitivePatch(expandedFacesTrg, neighbPatch().expandedPoints_));

		if (debug > 1)
		{
			Info << "Writing expanded geom patch as VTK" << endl;

			const Time& t = boundaryMesh().mesh().time();
			OFstream osSrc(t.path()/name() + "_extendedPatch.obj");
			OFstream osTrg(t.path()/neighbPatch().name() + "_extendedPatch.obj");
			meshTools::writeOBJ(osSrc, expandedFacesSrc, expandedPoints_);
			meshTools::writeOBJ(osTrg, expandedFacesTrg, neighbPatch().expandedPoints_);
		}

	    AMIPtr_->setRequireMatch(true);

	    // Construct/apply AMI interpolation to determine addressing and weights
	    AMIPtr_->upToDate() = false;
	    AMIPtr_->calculate(expandedPatchPtr_(), neighbPatch().expandedPatchPtr_());


        // Print some statistics
        const label nFace = returnReduce(size(), sumOp<label>());

        if (nFace)
        {
            scalarField srcWghtSum(size(), Zero);
            forAll(srcWghtSum, faceI)
            {
                srcWghtSum[faceI] = sum(AMIPtr_->srcWeights()[faceI]);
            }
            scalarField tgtWghtSum(neighbPatch().size(), Zero);
            forAll(tgtWghtSum, faceI)
            {
                tgtWghtSum[faceI] = sum(AMIPtr_->tgtWeights()[faceI]);
            }

            Info<< indent
                << "AMI: Patch " << name()
                << " sum(weights)"
                << " min:" << gMin(srcWghtSum)
                << " max:" << gMax(srcWghtSum)
                << " average:" << gAverage(srcWghtSum) << nl;
            Info<< indent
                << "AMI: Patch " << neighbPatch().name()
                << " sum(weights)"
                << " min:" << gMin(tgtWghtSum)
                << " max:" << gMax(tgtWghtSum)
                << " average:" << gAverage(tgtWghtSum) << nl;
        }

	    if (debug)
	    {
	        AMIPtr_->checkSymmetricWeights(true);
	    }
    }
}


void Foam::overlapAMIPolyPatch::calcTransforms()
{}


void Foam::overlapAMIPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    DebugInFunction << endl;

    // Flag AMI as needing update
    AMIPtr_->upToDate() = false;

    polyPatch::initGeometry(pBufs);

    calcTransforms();
}


void Foam::overlapAMIPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    DebugInFunction << endl;
}


void Foam::overlapAMIPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    DebugInFunction << endl;

    // See below. Clear out any local geometry
    primitivePatch::movePoints(p);

    // Note: processorPolyPatch::initMovePoints calls
    // processorPolyPatch::initGeometry which will trigger calculation of
    // patch faceCentres() and cell volumes...

    AMIPtr_->upToDate() = false;

    // Early calculation of transforms. See above.
    calcTransforms();
}


void Foam::overlapAMIPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    DebugInFunction << endl;

//    Note: not calling movePoints since this will undo our manipulations!
//    polyPatch::movePoints(pBufs, p);
/*
    polyPatch::movePoints
     -> primitivePatch::movePoints
        -> primitivePatch::clearGeom:
    deleteDemandDrivenData(localPointsPtr_);
    deleteDemandDrivenData(faceCentresPtr_);
    deleteDemandDrivenData(faceAreasPtr_);
    deleteDemandDrivenData(magFaceAreasPtr_);
    deleteDemandDrivenData(faceNormalsPtr_);
    deleteDemandDrivenData(pointNormalsPtr_);
*/
}


void Foam::overlapAMIPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    DebugInFunction << endl;

    polyPatch::initUpdateMesh(pBufs);
}


void Foam::overlapAMIPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    DebugInFunction << endl;

    // Note: this clears out cellCentres(), faceCentres() and faceAreas()
    polyPatch::updateMesh(pBufs);
}


void Foam::overlapAMIPolyPatch::clearGeom()
{
    DebugInFunction << endl;

    polyPatch::clearGeom();
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::overlapAMIPolyPatch::overlapAMIPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const transformType transform,
    const word& defaultAMIMethod
)
:
    coupledPolyPatch(name, size, start, index, bm, patchType, transform),
    rotationAxis_(Zero),
    nCopies_(0),
	expandedPatchPtr_(nullptr),
	expandedFaces_(size),
	expandedPoints_(size),
    nbrPatchName_(word::null),
    nbrPatchID_(-1),
    AMIPtr_(AMIInterpolation::New(defaultAMIMethod)),
    surfDict_(fileName("surface")),
    surfPtr_(nullptr)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::overlapAMIPolyPatch::overlapAMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const word& defaultAMIMethod
)
:
    coupledPolyPatch(name, dict, index, bm, patchType),
	rotationAxis_(dict.get<vector>("rotationAxis")),
	nCopies_(dict.get<label>("nCopies")),
	expandedPatchPtr_(nullptr),
	expandedFaces_(nCopies_*localFaces().size()),
	expandedPoints_(nCopies_*localPoints().size()),
    nbrPatchName_(dict.getOrDefault<word>("neighbourPatch", word::null)),
    coupleGroup_(dict),
    nbrPatchID_(-1),
    AMIPtr_
    (
        AMIInterpolation::New
        (
            dict.getOrDefault<word>("AMIMethod", defaultAMIMethod),
            dict,
            dict.getOrDefault("flipNormals", false)
        )
    ),
    surfDict_(dict.subOrEmptyDict("surface")),
    surfPtr_(nullptr)
{
    if (nbrPatchName_ == word::null && !coupleGroup_.valid())
    {
        FatalIOErrorInFunction(dict)
            << "No \"neighbourPatch\" or \"coupleGroup\" provided."
            << exit(FatalIOError);
    }

    if (nbrPatchName_ == name)
    {
        FatalIOErrorInFunction(dict)
            << "Neighbour patch name " << nbrPatchName_
            << " cannot be the same as this patch " << name
            << exit(FatalIOError);
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::overlapAMIPolyPatch::overlapAMIPolyPatch
(
    const overlapAMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
	rotationAxis_(pp.rotationAxis_),
	nCopies_(pp.nCopies_),
	expandedPatchPtr_(nullptr),
	expandedFaces_(pp.expandedFaces_),
	expandedPoints_(pp.expandedPoints_),
    nbrPatchName_(pp.nbrPatchName_),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    AMIPtr_(pp.AMIPtr_->clone()),
    surfDict_(pp.surfDict_),
    surfPtr_(nullptr)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::overlapAMIPolyPatch::overlapAMIPolyPatch
(
    const overlapAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
	rotationAxis_(pp.rotationAxis_),
	nCopies_(pp.nCopies_),
	expandedPatchPtr_(nullptr),
	expandedFaces_(pp.expandedFaces_),
	expandedPoints_(pp.expandedPoints_),
    nbrPatchName_(nbrPatchName),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    AMIPtr_(pp.AMIPtr_->clone()),
    surfDict_(pp.surfDict_),
    surfPtr_(nullptr)
{
    if (nbrPatchName_ == name())
    {
        FatalErrorInFunction
            << "Neighbour patch name " << nbrPatchName_
            << " cannot be the same as this patch " << name()
            << exit(FatalError);
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::overlapAMIPolyPatch::overlapAMIPolyPatch
(
    const overlapAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, mapAddressing, newStart),
	rotationAxis_(pp.rotationAxis_),
	nCopies_(pp.nCopies_),
	expandedPatchPtr_(nullptr),
	expandedFaces_(pp.expandedFaces_),
	expandedPoints_(pp.expandedPoints_),
    nbrPatchName_(pp.nbrPatchName_),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    AMIPtr_(pp.AMIPtr_->clone()),
    surfDict_(pp.surfDict_),
    surfPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::overlapAMIPolyPatch::neighbPatchID() const
{
    if (nbrPatchID_ == -1)
    {
        nbrPatchID_ = this->boundaryMesh().findPatchID(neighbPatchName());

        if (nbrPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal neighbourPatch name " << neighbPatchName()
                << nl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        // Check that it is a phaseLag AMI patch
        const overlapAMIPolyPatch& nbrPatch =
            refCast<const overlapAMIPolyPatch>
            (
                this->boundaryMesh()[nbrPatchID_]
            );

        if (nbrPatch.neighbPatchName() != name())
        {
            WarningInFunction
                << "Patch " << name()
                << " specifies neighbour patch " << neighbPatchName()
                << nl << " but that in return specifies "
                << nbrPatch.neighbPatchName() << endl;
        }
    }

    return nbrPatchID_;
}


bool Foam::overlapAMIPolyPatch::owner() const
{
    return index() < neighbPatchID();
}


const Foam::overlapAMIPolyPatch& Foam::overlapAMIPolyPatch::neighbPatch() const
{
    const polyPatch& pp = this->boundaryMesh()[neighbPatchID()];
    return refCast<const overlapAMIPolyPatch>(pp);
}


const Foam::AMIPatchToPatchInterpolation& Foam::overlapAMIPolyPatch::AMI() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "AMI interpolator only available to owner patch"
            << abort(FatalError);
    }

    if (!AMIPtr_->upToDate())
    {
        resetAMI();
    }

    return *AMIPtr_;
}


const Foam::primitivePatch& Foam::overlapAMIPolyPatch::expandedPatch() const
{
    if (!expandedPatchPtr_)
    {
        resetAMI();
    }

    return *expandedPatchPtr_;
}


bool Foam::overlapAMIPolyPatch::applyLowWeightCorrection() const
{
    if (owner())
    {
        return AMI().applyLowWeightCorrection();
    }
    else
    {
        return neighbPatch().AMI().applyLowWeightCorrection();
    }
}


void Foam::overlapAMIPolyPatch::transformPosition(pointField& l) const
{}


void Foam::overlapAMIPolyPatch::transformPosition
(
    point& l,
    const label facei
) const
{}


void Foam::overlapAMIPolyPatch::reverseTransformPosition
(
    point& l,
    const label facei
) const
{}


void Foam::overlapAMIPolyPatch::reverseTransformDirection
(
    vector& d,
    const label facei
) const
{
    if (!parallel())
    {
        const tensor& T =
        (
            reverseT().size() == 1
          ? reverseT()[0]
          : reverseT()[facei]
        );

        d = Foam::transform(T, d);
    }
}


void Foam::overlapAMIPolyPatch::calcGeometry
(
    const primitivePatch& referPatch,
    const pointField& thisCtrs,
    const vectorField& thisAreas,
    const pointField& thisCc,
    const pointField& nbrCtrs,
    const vectorField& nbrAreas,
    const pointField& nbrCc
)
{}


void Foam::overlapAMIPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{}


bool Foam::overlapAMIPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    return false;
}


Foam::label Foam::overlapAMIPolyPatch::pointFace
(
    const label facei,
    const vector& n,
    point& p
) const
{
    point prt(p);
    reverseTransformPosition(prt, facei);

    vector nrt(n);
    reverseTransformDirection(nrt, facei);

    label nbrFacei = -1;

    if (owner())
    {
        nbrFacei = AMI().tgtPointFace
        (
            expandedPatchPtr_(),
            neighbPatch().expandedPatchPtr_(),
            nrt,
            facei,
            prt
        );
    }
    else
    {
        nbrFacei = neighbPatch().AMI().srcPointFace
        (
            neighbPatch().expandedPatchPtr_(),
            expandedPatchPtr_(),
            nrt,
            facei,
            prt
        );
    }

    if (nbrFacei >= 0)
    {
        p = prt;
    }

    return nbrFacei;
}


Foam::tensor Foam::overlapAMIPolyPatch::RodriguesRotation
(
    const vector& rotationAxis,
    const scalar& rotationAngle,
    const bool inDegrees
) const
{
    tensor rotTensor;
    scalar theta = rotationAngle;

    if (inDegrees)
    {
        theta *= constant::mathematical::pi/180.0;
    }

    scalar sinTheta = sin(theta);
    scalar cosTheta = cos(theta);
    scalar oneMinusCosTheta = 1.0 - cosTheta;

    scalar magRotAxis = mag(rotationAxis);

    if (magRotAxis < SMALL)
    {
        FatalErrorIn
        (
            "tensor RodriguesRotation\n"
            "(\n"
            "    const vector& rotationAxis,\n"
            "    const scalar& rotationAngle\n"
            ")"
        )   << "Incorrectly defined axis: " << rotationAxis
            << abort(FatalError);
    }

    vector unitVector = rotationAxis/magRotAxis;

    scalar wx = unitVector.x();
    scalar wy = unitVector.y();
    scalar wz = unitVector.z();

    rotTensor.xx() = cosTheta + sqr(wx)*oneMinusCosTheta;
    rotTensor.yy() = cosTheta + sqr(wy)*oneMinusCosTheta;
    rotTensor.zz() = cosTheta + sqr(wz)*oneMinusCosTheta;

    rotTensor.xy() = wx*wy*oneMinusCosTheta - wz*sinTheta;
    rotTensor.xz() = wy*sinTheta + wx*wz*oneMinusCosTheta;

    rotTensor.yx() =  wz*sinTheta + wx*wy*oneMinusCosTheta;
    rotTensor.yz() = -wx*sinTheta + wy*wz*oneMinusCosTheta;

    rotTensor.zx() = -wy*sinTheta + wx*wz*oneMinusCosTheta;
    rotTensor.zy() =  wx*sinTheta + wy*wz*oneMinusCosTheta;

    return rotTensor;
}


void Foam::overlapAMIPolyPatch::write(Ostream& os) const
{
    coupledPolyPatch::write(os);
    if (!nbrPatchName_.empty())
    {
        os.writeEntry("neighbourPatch", nbrPatchName_);
    }
    coupleGroup_.write(os);

    AMIPtr_->write(os);

    if (!surfDict_.empty())
    {
        surfDict_.writeEntry(surfDict_.dictName(), os);
    }

    os.writeEntry("rotationAxis", rotationAxis_);
    os.writeEntry("nCopies", nCopies_);

}


// ************************************************************************* //
