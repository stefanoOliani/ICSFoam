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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline bool Foam::overlapAMIPolyPatch::canResetAMI() const
{
    return Pstream::parRun() || !boundaryMesh().mesh().time().processorCase();
}


inline const Foam::word& Foam::overlapAMIPolyPatch::neighbPatchName() const
{
    if (nbrPatchName_.empty())
    {
        // Try and use patchGroup to find samplePatch and sampleRegion
        label patchID = coupleGroup_.findOtherPatchID(*this);

        nbrPatchName_ = boundaryMesh()[patchID].name();
    }
    return nbrPatchName_;
}


inline const Foam::vector& Foam::overlapAMIPolyPatch::rotationAxis() const
{
    return rotationAxis_;
}


inline Foam::scalar Foam::overlapAMIPolyPatch::angle() const
{
    return 360.0/scalar(nCopies());
}


inline const Foam::label& Foam::overlapAMIPolyPatch::nCopies() const
{
    return nCopies_;
}


inline const Foam::scalarListList& Foam::overlapAMIPolyPatch::weights() const
{
    if (owner())
    {
        return AMI().srcWeights();
    }

    return neighbPatch().AMI().tgtWeights();
}


inline const Foam::scalarField& Foam::overlapAMIPolyPatch::weightsSum() const
{
    if (owner())
    {
        return AMI().srcWeightsSum();
    }

    return neighbPatch().AMI().tgtWeightsSum();
}


// ************************************************************************* //
