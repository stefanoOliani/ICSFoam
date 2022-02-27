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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline bool Foam::phaseLagAMIPolyPatch::canResetAMI() const
{
    return Pstream::parRun() || !boundaryMesh().mesh().time().processorCase();
}


inline const Foam::word& Foam::phaseLagAMIPolyPatch::neighbPatchName() const
{
    if (nbrPatchName_.empty())
    {
        // Try and use patchGroup to find samplePatch and sampleRegion
        label patchID = coupleGroup_.findOtherPatchID(*this);

        nbrPatchName_ = boundaryMesh()[patchID].name();
    }
    return nbrPatchName_;
}


inline Foam::label& Foam::phaseLagAMIPolyPatch::nTransformsFwd() const
{
    return nTransformsFwd_;
}


inline Foam::label& Foam::phaseLagAMIPolyPatch::nTransformsBwd() const
{
    return nTransformsBwd_;
}


inline Foam::label& Foam::phaseLagAMIPolyPatch::nTransforms() const
{
    return nTransforms_;
}


inline const Foam::scalarListList& Foam::phaseLagAMIPolyPatch::weights() const
{
    if (owner())
    {
        return AMI().srcWeights();
    }

    return neighbPatch().AMI().tgtWeights();
}


inline const Foam::scalarField& Foam::phaseLagAMIPolyPatch::weightsSum() const
{
    if (owner())
    {
        return AMI().srcWeightsSum();
    }

    return neighbPatch().AMI().tgtWeightsSum();
}


// ************************************************************************* //
