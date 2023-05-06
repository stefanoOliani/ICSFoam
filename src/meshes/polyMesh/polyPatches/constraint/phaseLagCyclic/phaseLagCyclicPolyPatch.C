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

#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "demandDrivenData.H"
#include "OFstream.H"
#include "matchPoints.H"
#include "edgeHashes.H"
#include "phaseLagCyclicPolyPatch.H"
#include "Time.H"
#include "transformField.H"
#include "SubField.H"
#include "unitConversion.H"

#include "IOHBZoneList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseLagCyclicPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, phaseLagCyclicPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, phaseLagCyclicPolyPatch, dictionary);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::phaseLagCyclicPolyPatch::phaseLagCyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const transformType transform
)
:
    cyclicPolyPatch(name, size, start, index, bm, patchType, transform),
	IBPA_(0),
	HBZoneName_(word::null),
	cylCoords_(false)
{}


Foam::phaseLagCyclicPolyPatch::phaseLagCyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicPolyPatch(name, dict, index, bm, patchType),
	IBPA_(dict.getOrDefault("IBPA", 0.0)),
	HBZoneName_(dict.getOrDefault("HBZoneName", word::null)),
	cylCoords_(dict.get<bool>("cylCoords"))
{}


Foam::phaseLagCyclicPolyPatch::phaseLagCyclicPolyPatch
(
    const phaseLagCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    cyclicPolyPatch(pp, bm),
	IBPA_(pp.IBPA_),
	HBZoneName_(pp.HBZoneName_),
	cylCoords_(pp.cylCoords_)
{}


Foam::phaseLagCyclicPolyPatch::phaseLagCyclicPolyPatch
(
    const phaseLagCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    cyclicPolyPatch(pp, bm, index, mapAddressing, newStart),
	IBPA_(pp.IBPA_),
	HBZoneName_(pp.HBZoneName_),
	cylCoords_(pp.cylCoords_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseLagCyclicPolyPatch::~phaseLagCyclicPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::phaseLagCyclicPolyPatch::timeLag() const
{
	const objectRegistry& allSubLevels = this->boundaryMesh().mesh().parent();
  	const objectRegistry& subLevel0 = allSubLevels.lookupObject<objectRegistry>("subTimeLevel0");
	const HBZoneList& HB = subLevel0.lookupObject<IOHBZoneList>("HBProperties");

	scalar IBPA = 0.0;
	label HBZoneInstance = -1;

	if (this->owner())
	{
		if (IBPA_ > 0)
		{
			IBPA = IBPA_;
		}
		else
		{
			IBPA = 2*constant::mathematical::pi + IBPA_;
		}

		forAll (HB, i)
		{
			if (HB[i].name() == HBZoneName_)
			{
				HBZoneInstance = i;
			}
		}
	}
	else
	{
		if (neighbPatch().IBPA() > 0)
		{
			IBPA = 2*constant::mathematical::pi - neighbPatch().IBPA();
		}
		else
		{
			IBPA = -neighbPatch().IBPA();
		}


		forAll (HB, i)
		{
			if (HB[i].name() == neighbPatch().HBZoneName())
			{
				HBZoneInstance = i;
			}
		}
	}

	if (IBPA < 0)
	{
	    FatalErrorInFunction
	        << "Negative time lag obtained."
	        << abort(FatalError);
	}

	const scalar& omega = HB[HBZoneInstance].omegaList()[1];

	return IBPA/omega;

}


void Foam::phaseLagCyclicPolyPatch::write(Ostream& os) const
{
    cyclicPolyPatch::write(os);

    os.writeEntryIfDifferent("IBPA", 0.0, IBPA_);
    os.writeEntryIfDifferent("HBZoneName", word::null, HBZoneName_);
    os.writeEntry("cylCoords", cylCoords_);

}

// ************************************************************************* //
