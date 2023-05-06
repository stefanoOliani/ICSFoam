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

#include "phaseLagAMIPolyPatch.H"
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

#include "IOHBZoneList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseLagAMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, phaseLagAMIPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, phaseLagAMIPolyPatch, dictionary);
}

const Foam::scalar Foam::phaseLagAMIPolyPatch::tolerance_ = 1e-10;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::phaseLagAMIPolyPatch::syncTransforms() const
{

	// Get the periodic patch
	const coupledPolyPatch& periodicPatch
	(
		refCast<const coupledPolyPatch>
		(
			boundaryMesh()[periodicPatchID()]
		)
	);

	// If there are any zero-sized periodic patches
	if (returnReduce((size() && !periodicPatch.size()), orOp<bool>()))
	{
		if (periodicPatch.separation().size() > 1)
		{
			FatalErrorInFunction
				<< "Periodic patch " << periodicPatchName_
				<< " has non-uniform separation vector "
				<< periodicPatch.separation()
				<< "This is not allowed inside " << type()
				<< " patch " << name()
				<< exit(FatalError);
		}

		if (periodicPatch.forwardT().size() > 1)
		{
			FatalErrorInFunction
				<< "Periodic patch " << periodicPatchName_
				<< " has non-uniform transformation tensor "
				<< periodicPatch.forwardT()
				<< "This is not allowed inside " << type()
				<< " patch " << name()
				<< exit(FatalError);
		}

		bool isParallel =
		(
			periodicPatch.size()
		 && periodicPatch.parallel()
		);
		reduce(isParallel, orOp<bool>());

		if (isParallel)
		{
			// Sync a list of separation vectors
			List<vectorField> sep(Pstream::nProcs());
			sep[Pstream::myProcNo()] = periodicPatch.separation();
			Pstream::gatherList(sep);
			Pstream::scatterList(sep);

			List<boolList> coll(Pstream::nProcs());
			coll[Pstream::myProcNo()] = periodicPatch.collocated();
			Pstream::gatherList(coll);
			Pstream::scatterList(coll);

			// If locally we have zero faces pick the first one that has a
			// separation vector
			if (!periodicPatch.size())
			{
				forAll(sep, procI)
				{
					if (sep[procI].size())
					{
						const_cast<vectorField&>
						(
							periodicPatch.separation()
						) = sep[procI];
						const_cast<boolList&>
						(
							periodicPatch.collocated()
						) = coll[procI];

						break;
					}
				}
			}
		}
		else
		{
			// Sync a list of forward and reverse transforms
			List<tensorField> forwardT(Pstream::nProcs());
			forwardT[Pstream::myProcNo()] = periodicPatch.forwardT();
			Pstream::gatherList(forwardT);
			Pstream::scatterList(forwardT);

			List<tensorField> reverseT(Pstream::nProcs());
			reverseT[Pstream::myProcNo()] = periodicPatch.reverseT();
			Pstream::gatherList(reverseT);
			Pstream::scatterList(reverseT);

			// If locally we have zero faces pick the first one that has a
			// transformation vector
			if (!periodicPatch.size())
			{
				forAll(forwardT, procI)
				{
					if (forwardT[procI].size())
					{
						const_cast<tensorField&>
						(
							periodicPatch.forwardT()
						) = forwardT[procI];
						const_cast<tensorField&>
						(
							periodicPatch.reverseT()
						) = reverseT[procI];

						break;
					}
				}
			}
		}
    }
}


void Foam::phaseLagAMIPolyPatch::writeOBJ
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
Foam::phaseLagAMIPolyPatch::surfPtr() const
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


void Foam::phaseLagAMIPolyPatch::resetAMI() const
{
    DebugInFunction << endl;

	const objectRegistry& allSubLevels = this->boundaryMesh().mesh().objectRegistry::parent();
  	const objectRegistry& subLevel0 = allSubLevels.lookupObject<objectRegistry>("subTimeLevel0");
	const HBZoneList& HB = subLevel0.lookupObject<IOHBZoneList>("HBProperties");
	const label nT = HB.selectedSnapshots().size();

    for (label i=0; i<nT; i++)
    {
    	word itrName = Foam::name(i);
    	word timeLevel = "subTimeLevel" + itrName;

    	if (timeLevel == this->boundaryMesh().mesh().name())
    	{
    		subTimeLevel_ = i;
    		neighbPatch().subTimeLevel_ = i;
    		break;
    	}
    }

    if (owner())
    {
        // Get the periodic patch
        const coupledPolyPatch& periodicPatchOwn
        (
            refCast<const coupledPolyPatch>
            (
                boundaryMesh()[periodicPatchID()]
            )
        );

        const coupledPolyPatch& periodicPatchNeigh
        (
            refCast<const coupledPolyPatch>
            (
                boundaryMesh()[neighbPatch().periodicPatchID()]
            )
        );

        // Synchronise the transforms
        syncTransforms();

        // Create copies of both patches' points
        pointField thisPoints0(localPoints());
        pointField nbrPoints0(neighbPatch().localPoints());

		pointField expandedPointsSrc(thisPoints0);
		faceList expandedFacesListSrc(localFaces());

		pointField expandedPointsTgt(nbrPoints0);
		faceList expandedFacesListTgt(neighbPatch().localFaces());

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

        // Create patches for all the points

        // Source patch at initial location
        primitivePatch thisPatch0
        (
            SubList<face>(localFaces(), size()),
            thisPoints0
        );
        // Target patch at initial location
        primitivePatch nbrPatch0
        (
            SubList<face>(neighbPatch().localFaces(), neighbPatch().size()),
            nbrPoints0
        );

        // Create another copy
        pointField thisPoints(thisPoints0);
        pointField nbrPoints(nbrPoints0);

        // Source patch that gets moved
        primitivePatch thisPatch
        (
            SubList<face>(localFaces(), size()),
            thisPoints
        );

        // Target patch that gets moved
        primitivePatch nbrPatch
        (
            SubList<face>(neighbPatch().localFaces(), neighbPatch().size()),
            nbrPoints
        );


        //TGT step

        label copyOffsetGeomPreSrc = 0;

        for (label i = 0; i < nTransformsBwd_; ++ i)
        {
        	faceList geomLocalFacesSrcPre(localFaces());

        	revTransformPositionPatch(thisPoints, periodicPatchOwn);

			expandedPointsSrc.append(thisPoints);

			thisPatch.movePoints(thisPoints);

			copyOffsetGeomPreSrc += thisPoints0.size();

			forAll (geomLocalFacesSrcPre, faceI)
			{
				face& curGeomFace = geomLocalFacesSrcPre[faceI];

				forAll (curGeomFace, fpI)
				{
					curGeomFace[fpI] += copyOffsetGeomPreSrc;
				}
			}

			expandedFacesListSrc.append(geomLocalFacesSrcPre);
        }

        pointField lastPointsOppDirSrc(thisPoints);
        thisPoints = thisPoints0;

        // Apply the stored number of periodic transforms
        for (label i = 0; i < nTransformsFwd_; ++ i)
        {
        	faceList geomLocalFacesSrcPre(localFaces());

            periodicPatchOwn.transformPosition(thisPoints);

            expandedPointsSrc.append(thisPoints);

            thisPatch.movePoints(thisPoints);

            copyOffsetGeomPreSrc += thisPoints0.size();

			forAll (geomLocalFacesSrcPre, faceI)
			{
				face& curGeomFace = geomLocalFacesSrcPre[faceI];

				forAll (curGeomFace, fpI)
				{
					curGeomFace[fpI] += copyOffsetGeomPreSrc;
				}
			}

			expandedFacesListSrc.append(geomLocalFacesSrcPre);
        }




        //SRC step

        label copyOffsetGeomPreTgt = 0;

        for (label i = 0; i < neighbPatch().nTransformsBwd(); ++ i)
        {
        	faceList geomLocalFacesTgtPre(neighbPatch().localFaces());

        	neighbPatch().revTransformPositionPatch(nbrPoints, periodicPatchNeigh);

			expandedPointsTgt.append(nbrPoints);

			nbrPatch.movePoints(nbrPoints);

			copyOffsetGeomPreTgt += nbrPoints0.size();

			forAll (geomLocalFacesTgtPre, faceI)
			{
				face& curGeomFace = geomLocalFacesTgtPre[faceI];

				forAll (curGeomFace, fpI)
				{
					curGeomFace[fpI] += copyOffsetGeomPreTgt;
				}
			}

			expandedFacesListTgt.append(geomLocalFacesTgtPre);
        }

        pointField lastPointsOppDirTrg(nbrPoints);
        nbrPoints = nbrPoints0;

        // Apply the stored number of periodic transforms
        for (label i = 0; i < neighbPatch().nTransformsFwd(); ++ i)
        {
        	faceList geomLocalFacesTgtPre(neighbPatch().localFaces());

        	periodicPatchNeigh.transformPosition(nbrPoints);

			expandedPointsTgt.append(nbrPoints);

			nbrPatch.movePoints(nbrPoints);

			copyOffsetGeomPreTgt += nbrPoints0.size();

			forAll (geomLocalFacesTgtPre, faceI)
			{
				face& curGeomFace = geomLocalFacesTgtPre[faceI];

				forAll (curGeomFace, fpI)
				{
					curGeomFace[fpI] += copyOffsetGeomPreTgt;
				}
			}

			expandedFacesListTgt.append(geomLocalFacesTgtPre);
        }

        pointField	expandedPointsSrcPre(expandedPointsSrc);
        pointField	expandedPointsTrgPre(expandedPointsTgt);

        faceList expandedFacesListSrcPre(expandedFacesListSrc);
        faceList expandedFacesListTgtPre(expandedFacesListTgt);

		SubList<face> expandedFacesSrcPre(expandedFacesListSrcPre);
		primitivePatch expandedPatchSrcPre(expandedFacesSrcPre, expandedPointsSrcPre);

		SubList<face> expandedFacesTrgPre(expandedFacesListTgtPre);
		primitivePatch expandedPatchTrgPre(expandedFacesTrgPre, expandedPointsTrgPre);

        // Construct a new AMI interpolation between the initial patch locations
        AMIPtr_->setRequireMatch(false);
	    AMIPtr_->upToDate() = false;
	    AMIPtr_->calculate(expandedPatchSrcPre, nbrPatch0);

        // Number of geometry replications
        label iter(0);

        if (ownStr.valid())
        {
            writeOBJ(thisPatch0, ownStr());
        }
        if (neiStr.valid())
        {
            writeOBJ(nbrPatch0, neiStr());
        }

        // Weight sum averages
        scalar srcSum(gAverage(AMIPtr_->srcWeightsSum()));
        scalar tgtSum(gAverage(AMIPtr_->tgtWeightsSum()));

        // Direction (or rather side of AMI : this or nbr patch) of
        // geometry replication: start always in positive direction
        bool directionSrc = 1;
        bool directionTgt = 1;

        // Increase in the source and target weight sum for the last iteration in the
        // opposite direction. If the current increase is less than this, the
        // direction (= side of AMI to transform) is reversed.

        scalar srcSumDiff = 0;
        scalar trgSumDiff = 0;

        DebugInFunction
            << "patch:" << name()
            << " srcSum:" << srcSum
            << " tgtSum:" << tgtSum
            << endl;

        // TGT step, transform the source until the target is completely covered


        while
        (
            (iter < maxIter_)
         && (
                gMin(AMIPtr_->tgtWeightsSum()) < matchTolerance()
            )
        )
        {
    		faceList geomLocalFacesSrc(localFaces());

            if (directionSrc)
            {
                periodicPatchOwn.transformPosition(thisPoints);

                expandedPointsSrc.append(thisPoints);

                thisPatch.movePoints(thisPoints);

            	const label copyOffsetGeom = (nTransforms_+1)*thisPoints0.size();

    			forAll (geomLocalFacesSrc, faceI)
    			{
    				face& curGeomFace = geomLocalFacesSrc[faceI];

    				forAll (curGeomFace, fpI)
    				{
    					curGeomFace[fpI] += copyOffsetGeom;
    				}
    			}

    			expandedFacesListSrc.append(geomLocalFacesSrc);

    			SubList<face> expandedFacesSrc(expandedFacesListSrc);
    			primitivePatch expandedPatchSrc(expandedFacesSrc, expandedPointsSrc);

    			AMIPtr_->upToDate() = false;
    			AMIPtr_->setRequireMatch(false);
                AMIPtr_->calculate(expandedPatchSrc, nbrPatch0);

                nTransformsFwd_ ++;

                if (ownStr.valid())
                {
                    writeOBJ(thisPatch, ownStr());
                }
            }
            else
            {
            	revTransformPositionPatch(thisPoints, periodicPatchOwn);

                expandedPointsSrc.append(thisPoints);

                thisPatch.movePoints(thisPoints);

            	const label copyOffsetGeom = (nTransforms_+1)*thisPoints0.size();

    			forAll (geomLocalFacesSrc, faceI)
    			{
    				face& curGeomFace = geomLocalFacesSrc[faceI];

    				forAll (curGeomFace, fpI)
    				{
    					curGeomFace[fpI] += copyOffsetGeom;
    				}
    			}

    			expandedFacesListSrc.append(geomLocalFacesSrc);

    			SubList<face> expandedFacesSrc(expandedFacesListSrc);
    			primitivePatch expandedPatchSrc(expandedFacesSrc, expandedPointsSrc);

    			AMIPtr_->upToDate() = false;
    			AMIPtr_->setRequireMatch(false);
                AMIPtr_->calculate(expandedPatchSrc, nbrPatch0);

                nTransformsBwd_ ++;

                if (ownStr.valid())
                {
                    writeOBJ(thisPatch, ownStr());
                }
            }

            const scalar trgSumNew = gAverage(AMIPtr_->tgtWeightsSum());
            const scalar trgSumDiffNew = trgSumNew - tgtSum;

            if (trgSumDiffNew < trgSumDiff || trgSumDiffNew < SMALL)
            {
                directionSrc = !directionSrc;

                trgSumDiff = trgSumDiffNew;

                pointField thisPointsTmp = thisPoints;

                //If changing direction, restart from last
                //position in the opposite direction
                thisPoints = lastPointsOppDirSrc;

                lastPointsOppDirSrc = thisPointsTmp;
            }

            tgtSum = trgSumNew;

            nTransforms_ = nTransformsFwd_ + nTransformsBwd_;

            ++iter;

            DebugInFunction
                << "patch:" << name()
                << " iteration:" << iter
                << " srcSum:" << srcSum
                << " tgtSum:" << tgtSum
                << " directionSrc:" << directionSrc
                << endl;
        }

        AMIPtr_->upToDate() = false;
        AMIPtr_->setRequireMatch(false);
        AMIPtr_->calculate(thisPatch0, expandedPatchTrgPre);

        // SRC step, transform the target until the source is completely covered

        while
        (
            (iter < maxIter_)
         && (
        		 gMin(AMIPtr_->srcWeightsSum()) < matchTolerance()
            )
        )
        {
    		faceList geomLocalFacesTgt(neighbPatch().localFaces());

            if (directionTgt)
            {
            	periodicPatchNeigh.transformPosition(nbrPoints);

                expandedPointsTgt.append(nbrPoints);

                nbrPatch.movePoints(nbrPoints);

            	const label copyOffsetGeom = (neighbPatch().nTransforms() + 1)*nbrPoints0.size();

    			forAll (geomLocalFacesTgt, faceI)
    			{
    				face& curGeomFace = geomLocalFacesTgt[faceI];

    				forAll (curGeomFace, fpI)
    				{
    					curGeomFace[fpI] += copyOffsetGeom;
    				}
    			}

    			expandedFacesListTgt.append(geomLocalFacesTgt);

    			SubList<face> expandedFacesTgt(expandedFacesListTgt);
    			primitivePatch expandedPatchTgt(expandedFacesTgt, expandedPointsTgt);

    			AMIPtr_->upToDate() = false;
    			AMIPtr_->setRequireMatch(false);
                AMIPtr_->calculate(thisPatch0, expandedPatchTgt);

                neighbPatch().nTransformsFwd() ++;

                if (neiStr.valid())
                {
                    writeOBJ(nbrPatch, neiStr());
                }
            }
            else
            {
            	neighbPatch().revTransformPositionPatch(nbrPoints, periodicPatchNeigh);

                expandedPointsTgt.append(nbrPoints);

                thisPatch.movePoints(nbrPoints);

            	const label copyOffsetGeom = (neighbPatch().nTransforms() + 1)*nbrPoints0.size();

    			forAll (geomLocalFacesTgt, faceI)
    			{
    				face& curGeomFace = geomLocalFacesTgt[faceI];

    				forAll (curGeomFace, fpI)
    				{
    					curGeomFace[fpI] += copyOffsetGeom;
    				}
    			}

    			expandedFacesListTgt.append(geomLocalFacesTgt);

    			SubList<face> expandedFacesTgt(expandedFacesListTgt);
    			primitivePatch expandedPatchTgt(expandedFacesTgt, expandedPointsTgt);

    			AMIPtr_->upToDate() = false;
    			AMIPtr_->setRequireMatch(false);
                AMIPtr_->calculate(thisPatch0, expandedPatchTgt);

                neighbPatch().nTransformsBwd() ++;

                if (neiStr.valid())
                {
                    writeOBJ(nbrPatch, neiStr());
                }
            }

            const scalar srcSumNew = gAverage(AMIPtr_->srcWeightsSum());
            const scalar srcSumDiffNew = srcSumNew - srcSum;

            if (srcSumDiffNew < srcSumDiff || srcSumDiffNew < SMALL)
            {
                directionTgt = !directionTgt;

                srcSumDiff = srcSumDiffNew;

                pointField nbrPointsTmp = nbrPoints;

                //If changing direction, restart from last
                //position in the opposite direction
                nbrPoints = lastPointsOppDirTrg;

                lastPointsOppDirTrg = nbrPointsTmp;
            }

            srcSum = srcSumNew;

            neighbPatch().nTransforms() = neighbPatch().nTransformsFwd() + neighbPatch().nTransformsBwd();

            ++iter;

            DebugInFunction
                << "patch:" << neighbPatch().name()
                << " iteration:" << iter
                << " srcSum:" << srcSum
                << " tgtSum:" << tgtSum
                << " directionSrc:" << directionTgt
                << endl;
        }




        //Reorder the interpolation

        thisPoints = thisPoints0;
        nbrPoints = nbrPoints0;

        expandedPoints_ = thisPoints0;
        expandedFaces_ = localFaces();

        neighbPatch().expandedPoints_ = nbrPoints0;
        neighbPatch().expandedFaces_ = neighbPatch().localFaces();

        //TGT step
		label copyOffsetGeomPostSrc = 0;

		for (label i = 0; i < nTransformsBwd_; ++ i)
		{
			faceList geomLocalFacesSrcPost(localFaces());

			revTransformPositionPatch(thisPoints, periodicPatchOwn);

			expandedPoints_.append(thisPoints);

			thisPatch.movePoints(thisPoints);

			copyOffsetGeomPostSrc += thisPoints0.size();

			forAll (geomLocalFacesSrcPost, faceI)
			{
				face& curGeomFace = geomLocalFacesSrcPost[faceI];

				forAll (curGeomFace, fpI)
				{
					curGeomFace[fpI] += copyOffsetGeomPostSrc;
				}
			}

			expandedFaces_.append(geomLocalFacesSrcPost);
		}

		thisPoints = thisPoints0;

		// Apply the stored number of periodic transforms
		for (label i = 0; i < nTransformsFwd_; ++ i)
		{
			faceList geomLocalFacesSrcPost(localFaces());

			periodicPatchOwn.transformPosition(thisPoints);

			expandedPoints_.append(thisPoints);

			thisPatch.movePoints(thisPoints);

			copyOffsetGeomPostSrc += thisPoints0.size();

			forAll (geomLocalFacesSrcPost, faceI)
			{
				face& curGeomFace = geomLocalFacesSrcPost[faceI];

				forAll (curGeomFace, fpI)
				{
					curGeomFace[fpI] += copyOffsetGeomPostSrc;
				}
			}

			expandedFaces_.append(geomLocalFacesSrcPost);
		}

		//SRC step

		label copyOffsetGeomPostTgt = 0;

		for (label i = 0; i < neighbPatch().nTransformsBwd(); ++ i)
		{
			faceList geomLocalFacesTgtPost(neighbPatch().localFaces());

			neighbPatch().revTransformPositionPatch(nbrPoints, periodicPatchNeigh);

			neighbPatch().expandedPoints_.append(nbrPoints);

			nbrPatch.movePoints(nbrPoints);

			copyOffsetGeomPostTgt += nbrPoints0.size();

			forAll (geomLocalFacesTgtPost, faceI)
			{
				face& curGeomFace = geomLocalFacesTgtPost[faceI];

				forAll (curGeomFace, fpI)
				{
					curGeomFace[fpI] += copyOffsetGeomPostTgt;
				}
			}

			neighbPatch().expandedFaces_.append(geomLocalFacesTgtPost);
		}

		nbrPoints = nbrPoints0;

		// Apply the stored number of periodic transforms
		for (label i = 0; i < neighbPatch().nTransformsFwd(); ++ i)
		{
			faceList geomLocalFacesTgtPost(neighbPatch().localFaces());

			periodicPatchNeigh.transformPosition(nbrPoints);

			neighbPatch().expandedPoints_.append(nbrPoints);

			nbrPatch.movePoints(nbrPoints);

			copyOffsetGeomPostTgt += nbrPoints0.size();

			forAll (geomLocalFacesTgtPost, faceI)
			{
				face& curGeomFace = geomLocalFacesTgtPost[faceI];

				forAll (curGeomFace, fpI)
				{
					curGeomFace[fpI] += copyOffsetGeomPostTgt;
				}
			}

			neighbPatch().expandedFaces_.append(geomLocalFacesTgtPost);
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

	    AMIPtr_->setRequireMatch(false);

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


void Foam::phaseLagAMIPolyPatch::calcTransforms()
{}


void Foam::phaseLagAMIPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    DebugInFunction << endl;

    // Flag AMI as needing update
    AMIPtr_->upToDate() = false;

    polyPatch::initGeometry(pBufs);

    calcTransforms();
}


void Foam::phaseLagAMIPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    DebugInFunction << endl;
}


void Foam::phaseLagAMIPolyPatch::initMovePoints
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


void Foam::phaseLagAMIPolyPatch::movePoints
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


void Foam::phaseLagAMIPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    DebugInFunction << endl;

    polyPatch::initUpdateMesh(pBufs);
}


void Foam::phaseLagAMIPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    DebugInFunction << endl;

    // Note: this clears out cellCentres(), faceCentres() and faceAreas()
    polyPatch::updateMesh(pBufs);
}


void Foam::phaseLagAMIPolyPatch::clearGeom()
{
    DebugInFunction << endl;

    polyPatch::clearGeom();
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::phaseLagAMIPolyPatch::phaseLagAMIPolyPatch
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
    periodicPatchName_(word::null),
    periodicPatchID_(-1),
	HBZoneName_(),
	subTimeLevel_(0),
	IBPA_(),
	cylCoords_(false),
	nTransformsFwd_(0),
	nTransformsBwd_(0),
    nTransforms_(0),
	expandedPatchPtr_(nullptr),
	expandedFaces_(size),
	expandedPoints_(size),
    nSectors_(0),
    maxIter_(36),
    nbrPatchName_(word::null),
    nbrPatchID_(-1),
    AMIPtr_(AMIInterpolation::New(defaultAMIMethod)),
    surfDict_(fileName("surface")),
    surfPtr_(nullptr)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::phaseLagAMIPolyPatch::phaseLagAMIPolyPatch
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
    periodicPatchName_(dict.lookup("periodicPatch")),
    periodicPatchID_(-1),
	HBZoneName_(dict.lookup("HBZoneName")),
	subTimeLevel_(0),
	IBPA_(dict.getOrDefault<scalar>("IBPA", 0)),
	cylCoords_(dict.get<bool>("cylCoords")),
	nTransformsFwd_(0),
	nTransformsBwd_(0),
    nTransforms_(0),
	expandedPatchPtr_(nullptr),
	expandedFaces_(localFaces().size()),
	expandedPoints_(localPoints().size()),
    nSectors_(dict.getOrDefault<label>("nSectors", 0)),
    maxIter_(dict.getOrDefault<label>("maxIter", 36)),
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


Foam::phaseLagAMIPolyPatch::phaseLagAMIPolyPatch
(
    const phaseLagAMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    periodicPatchName_(pp.periodicPatchName_),
    periodicPatchID_(-1),
	HBZoneName_(pp.HBZoneName_),
	subTimeLevel_(pp.subTimeLevel_),
	IBPA_(pp.IBPA_),
	cylCoords_(pp.cylCoords_),
	nTransformsFwd_(0),
	nTransformsBwd_(0),
    nTransforms_(0),
	expandedPatchPtr_(nullptr),
	expandedFaces_(pp.expandedFaces_),
	expandedPoints_(pp.expandedPoints_),
    nSectors_(pp.nSectors_),
    maxIter_(pp.maxIter_),
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


Foam::phaseLagAMIPolyPatch::phaseLagAMIPolyPatch
(
    const phaseLagAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    periodicPatchName_(pp.periodicPatchName_),
    periodicPatchID_(-1),
	HBZoneName_(pp.HBZoneName_),
	subTimeLevel_(pp.subTimeLevel_),
	IBPA_(pp.IBPA_),
	cylCoords_(pp.cylCoords_),
	nTransformsFwd_(0),
	nTransformsBwd_(0),
    nTransforms_(0),
	expandedPatchPtr_(nullptr),
	expandedFaces_(pp.expandedFaces_),
	expandedPoints_(pp.expandedPoints_),
    nSectors_(pp.nSectors_),
    maxIter_(pp.maxIter_),
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


Foam::phaseLagAMIPolyPatch::phaseLagAMIPolyPatch
(
    const phaseLagAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, mapAddressing, newStart),
    periodicPatchName_(pp.periodicPatchName_),
    periodicPatchID_(-1),
	HBZoneName_(pp.HBZoneName_),
	subTimeLevel_(pp.subTimeLevel_),
	IBPA_(pp.IBPA_),
	cylCoords_(pp.cylCoords_),
	nTransformsFwd_(0),
	nTransformsBwd_(0),
    nTransforms_(0),
	expandedPatchPtr_(nullptr),
	expandedFaces_(pp.expandedFaces_),
	expandedPoints_(pp.expandedPoints_),
    nSectors_(pp.nSectors_),
    maxIter_(pp.maxIter_),
    nbrPatchName_(pp.nbrPatchName_),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    AMIPtr_(pp.AMIPtr_->clone()),
    surfDict_(pp.surfDict_),
    surfPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::phaseLagAMIPolyPatch::periodicPatchID() const
{
    if (periodicPatchName_ == word::null)
    {
        periodicPatchID_ = -1;

        return periodicPatchID_;
    }

    if (periodicPatchID_ == -1)
    {
        periodicPatchID_ = this->boundaryMesh().findPatchID(periodicPatchName_);

        if (periodicPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal periodicPatch name " << periodicPatchName_
                << nl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        // Check that it is a coupled patch
        refCast<const coupledPolyPatch>
        (
            this->boundaryMesh()[periodicPatchID_]
        );
    }

    return periodicPatchID_;
}


Foam::label Foam::phaseLagAMIPolyPatch::neighbPatchID() const
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
        const phaseLagAMIPolyPatch& nbrPatch =
            refCast<const phaseLagAMIPolyPatch>
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


bool Foam::phaseLagAMIPolyPatch::owner() const
{
    return index() < neighbPatchID();
}


const Foam::phaseLagAMIPolyPatch& Foam::phaseLagAMIPolyPatch::neighbPatch() const
{
    const polyPatch& pp = this->boundaryMesh()[neighbPatchID()];
    return refCast<const phaseLagAMIPolyPatch>(pp);
}


const Foam::AMIPatchToPatchInterpolation& Foam::phaseLagAMIPolyPatch::AMI() const
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


bool Foam::phaseLagAMIPolyPatch::applyLowWeightCorrection() const
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


Foam::scalar Foam::phaseLagAMIPolyPatch::timeLag(bool fwd) const
{
	const objectRegistry& allSubLevels = this->boundaryMesh().mesh().parent();
  	const objectRegistry& subLevel0 = allSubLevels.lookupObject<objectRegistry>("subTimeLevel0");
	const HBZoneList& HB = subLevel0.lookupObject<IOHBZoneList>("HBProperties");

	scalar IBPA = 0.0;
	label HBZoneInstance = -1;

	if (fwd)
	{
		if (IBPA_ > 0)
		{
			IBPA = IBPA_;
		}
		else
		{
			IBPA = 2*constant::mathematical::pi + IBPA_;
		}
	}
	else
	{
		if (IBPA_ > 0)
		{
			IBPA = 2*constant::mathematical::pi - IBPA_;
		}
		else
		{
			IBPA = - IBPA_;
		}
	}

	forAll (HB, i)
	{
		if (HB[i].name() == HBZoneName_)
		{
			HBZoneInstance = i;
		}
	}

	if (IBPA < 0)
	{
	    FatalErrorInFunction
	        << "Negative time lag obtained."
	        << abort(FatalError);
	}

	if (!HB[HBZoneInstance].active())
	{
		return 0.0;
	}

	const scalar& omega = HB[HBZoneInstance].omegaList()[1];

	return IBPA/omega;

}


void Foam::phaseLagAMIPolyPatch::transformPosition(pointField& l) const
{}


void Foam::phaseLagAMIPolyPatch::transformPosition
(
    point& l,
    const label facei
) const
{}


void Foam::phaseLagAMIPolyPatch::reverseTransformPosition
(
    point& l,
    const label facei
) const
{}


void Foam::phaseLagAMIPolyPatch::reverseTransformDirection
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


void Foam::phaseLagAMIPolyPatch::revTransformPositionPatch
(
	pointField& l,
	const coupledPolyPatch& perioPatch
) const
{
    if (!perioPatch.parallel())
    {
    	//Please notice that if the transform is rotational
    	//and the rotation center is different from (0 0 0)
    	//the formula below is not correct. Unfortunately
    	//rotationCenter_ is not defined at the
    	//coupled poly patch level
		l = Foam::transform(perioPatch.reverseT(), l);
    }
    else if (perioPatch.separated())
    {
        // transformPosition gets called on the receiving side,
        // separation gets calculated on the sending side so subtract
    	// here we add since it is the reverse transform

        const vectorField& s = perioPatch.separation();
        if (s.size() == 1)
        {
            forAll(l, i)
            {
                l[i] += s[0];
            }
        }
        else
        {
            l += s;
        }
    }
}


void Foam::phaseLagAMIPolyPatch::calcGeometry
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


void Foam::phaseLagAMIPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{}


bool Foam::phaseLagAMIPolyPatch::order
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


Foam::label Foam::phaseLagAMIPolyPatch::pointFace
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


void Foam::phaseLagAMIPolyPatch::write(Ostream& os) const
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

    os.writeEntry("periodicPatch", periodicPatchName_);
    os.writeEntry("HBZoneName", HBZoneName_);
    os.writeEntry("IBPA", IBPA_);
    os.writeEntry("cylCoords", cylCoords_);
    os.writeEntryIfDifferent<label>("nTransforms", 0, nTransforms_);
    os.writeEntryIfDifferent<label>("nSectors", 0, nSectors_);
    os.writeEntryIfDifferent<label>("maxIter", 36, maxIter_);
}


// ************************************************************************* //
