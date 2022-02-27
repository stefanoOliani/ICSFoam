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

#include "MRFTranslatingZone.H"

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"
#include "faceSet.H"
#include "geometricOneField.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MRFTranslatingZone, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::MRFTranslatingZone::setMRFFaces()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Type per face:
    //  0:not in zone
    //  1:moving with frame
    //  2:other
    labelList faceType(mesh_.nFaces(), Zero);

    // Determine faces in cell zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (without constructing cells)

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    // Cells in zone
    boolList zoneCell(mesh_.nCells(), false);

    if (cellZoneID_ != -1)
    {
        const labelList& cellLabels = mesh_.cellZones()[cellZoneID_];
        forAll(cellLabels, i)
        {
            zoneCell[cellLabels[i]] = true;
        }
    }


    label nZoneFaces = 0;

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (zoneCell[own[facei]] || zoneCell[nei[facei]])
        {
            faceType[facei] = 1;
            nZoneFaces++;
        }
    }


    labelHashSet excludedPatches(excludedPatchLabels_);

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled() || excludedPatches.found(patchi))
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;

                if (zoneCell[own[facei]])
                {
                    faceType[facei] = 2;
                    nZoneFaces++;
                }
            }
        }
        else if (!isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;

                if (zoneCell[own[facei]])
                {
                    faceType[facei] = 1;
                    nZoneFaces++;
                }
            }
        }
    }

    // Synchronize the faceType across processor patches
    syncTools::syncFaceList(mesh_, faceType, maxEqOp<label>());

    // Now we have for faceType:
    //  0   : face not in cellZone
    //  1   : internal face or normal patch face
    //  2   : coupled patch face or excluded patch face

    // Sort into lists per patch.

    internalFaces_.setSize(mesh_.nFaces());
    label nInternal = 0;

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (faceType[facei] == 1)
        {
            internalFaces_[nInternal++] = facei;
        }
    }
    internalFaces_.setSize(nInternal);

    labelList nIncludedFaces(patches.size(), Zero);
    labelList nExcludedFaces(patches.size(), Zero);

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        forAll(pp, patchFacei)
        {
            label facei = pp.start() + patchFacei;

            if (faceType[facei] == 1)
            {
                nIncludedFaces[patchi]++;
            }
            else if (faceType[facei] == 2)
            {
                nExcludedFaces[patchi]++;
            }
        }
    }

    includedFaces_.setSize(patches.size());
    excludedFaces_.setSize(patches.size());
    forAll(nIncludedFaces, patchi)
    {
        includedFaces_[patchi].setSize(nIncludedFaces[patchi]);
        excludedFaces_[patchi].setSize(nExcludedFaces[patchi]);
    }
    nIncludedFaces = 0;
    nExcludedFaces = 0;

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        forAll(pp, patchFacei)
        {
            label facei = pp.start() + patchFacei;

            if (faceType[facei] == 1)
            {
                includedFaces_[patchi][nIncludedFaces[patchi]++] = patchFacei;
            }
            else if (faceType[facei] == 2)
            {
                excludedFaces_[patchi][nExcludedFaces[patchi]++] = patchFacei;
            }
        }
    }


    if (debug)
    {
        faceSet internalFaces(mesh_, "internalFaces", internalFaces_);
        Pout<< "Writing " << internalFaces.size()
            << " internal faces in MRF zone to faceSet "
            << internalFaces.name() << endl;
        internalFaces.write();

        faceSet MRFFaces(mesh_, "includedFaces", 100);
        forAll(includedFaces_, patchi)
        {
            forAll(includedFaces_[patchi], i)
            {
                label patchFacei = includedFaces_[patchi][i];
                MRFFaces.insert(patches[patchi].start()+patchFacei);
            }
        }
        Pout<< "Writing " << MRFFaces.size()
            << " patch faces in MRF zone to faceSet "
            << MRFFaces.name() << endl;
        MRFFaces.write();

        faceSet excludedFaces(mesh_, "excludedFaces", 100);
        forAll(excludedFaces_, patchi)
        {
            forAll(excludedFaces_[patchi], i)
            {
                label patchFacei = excludedFaces_[patchi][i];
                excludedFaces.insert(patches[patchi].start()+patchFacei);
            }
        }
        Pout<< "Writing " << excludedFaces.size()
            << " faces in MRF zone with special handling to faceSet "
            << excludedFaces.name() << endl;
        excludedFaces.write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFTranslatingZone::MRFTranslatingZone
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    mesh_(mesh),
    name_(name),
    coeffs_(dict),
    active_(coeffs_.getOrDefault("active", true)),
    cellZoneName_(cellZoneName),
    cellZoneID_(),
    excludedPatchNames_
    (
        coeffs_.getOrDefault<wordRes>("nonMovingPatches", wordRes())
    ),
    frameVel_(Function1<vector>::New("frameVelocity", coeffs_))
{
    if (cellZoneName_ == word::null)
    {
        coeffs_.readEntry("cellZone", cellZoneName_);
    }

    if (!active_)
    {
        cellZoneID_ = -1;
    }
    else
    {
        cellZoneID_ = mesh_.cellZones().findZoneID(cellZoneName_);

        const labelHashSet excludedPatchSet
        (
            mesh_.boundaryMesh().patchSet(excludedPatchNames_)
        );

        excludedPatchLabels_.setSize(excludedPatchSet.size());

        label i = 0;
        for (const label patchi : excludedPatchSet)
        {
            excludedPatchLabels_[i++] = patchi;
        }

        bool cellZoneFound = (cellZoneID_ != -1);

        reduce(cellZoneFound, orOp<bool>());

        if (!cellZoneFound)
        {
            FatalErrorInFunction
                << "cannot find MRF cellZone " << cellZoneName_
                << exit(FatalError);
        }

        setMRFFaces();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::MRFTranslatingZone::frameVelocity() const
{
    return frameVel_->value(mesh_.time().timeOutputValue());
}


void Foam::MRFTranslatingZone::addFrameVelocity(volVectorField& UFrame) const
{
    const vector Velocity = frameVelocity();

    // Set omega in all cells of the rotating cell zone
    const labelList& cells = mesh_.cellZones()[cellZoneID_];

    forAll (cells, i)
    {
        UFrame[cells[i]] = Velocity;
    }
}


void Foam::MRFTranslatingZone::makeRelative(volVectorField& U) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const vector Velocity = this->frameVelocity();

    const labelList& cells = mesh_.cellZones()[cellZoneID_];

    forAll(cells, i)
    {
        label celli = cells[i];
        U[celli] -= Velocity;
    }

    // Included patches

    volVectorField::Boundary& Ubf = U.boundaryFieldRef();

    forAll(includedFaces_, patchi)
    {
        forAll(includedFaces_[patchi], i)
        {
            label patchFacei = includedFaces_[patchi][i];
            Ubf[patchi][patchFacei] = Zero;
        }
    }

    // Excluded patches
    forAll(excludedFaces_, patchi)
    {
        forAll(excludedFaces_[patchi], i)
        {
            label patchFacei = excludedFaces_[patchi][i];
            Ubf[patchi][patchFacei] -= Velocity;
        }
    }
}


void Foam::MRFTranslatingZone::makeRelative(surfaceScalarField& phi) const
{
    makeRelativeRhoFlux(geometricOneField(), phi);
}


void Foam::MRFTranslatingZone::makeRelative(FieldField<fvsPatchField, scalar>& phi) const
{
    makeRelativeRhoFlux(oneFieldField(), phi);
}


void Foam::MRFTranslatingZone::makeRelative(Field<scalar>& phi, const label patchi) const
{
    makeRelativeRhoFlux(oneField(), phi, patchi);
}


void Foam::MRFTranslatingZone::makeRelative
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    makeRelativeRhoFlux(rho, phi);
}


void Foam::MRFTranslatingZone::makeAbsolute(volVectorField& U) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const vector Velocity = this->frameVelocity();

    const labelList& cells = mesh_.cellZones()[cellZoneID_];

    forAll(cells, i)
    {
        label celli = cells[i];
        U[celli] += Velocity;
    }

    // Included patches
    volVectorField::Boundary& Ubf = U.boundaryFieldRef();

    forAll(includedFaces_, patchi)
    {
        forAll(includedFaces_[patchi], i)
        {
            label patchFacei = includedFaces_[patchi][i];
            Ubf[patchi][patchFacei] = Velocity;
        }
    }

    // Excluded patches
    forAll(excludedFaces_, patchi)
    {
        forAll(excludedFaces_[patchi], i)
        {
            label patchFacei = excludedFaces_[patchi][i];
            Ubf[patchi][patchFacei] += Velocity;
        }
    }
}


void Foam::MRFTranslatingZone::makeAbsolute(surfaceScalarField& phi) const
{
    makeAbsoluteRhoFlux(geometricOneField(), phi);
}


void Foam::MRFTranslatingZone::makeAbsolute
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    makeAbsoluteRhoFlux(rho, phi);
}


void Foam::MRFTranslatingZone::faceU
(
    surfaceVectorField& zoneFaceU
) const
{
    const vector& Velocity = this->frameVelocity();

    // Internal faces
    forAll(internalFaces_, i)
    {
        label facei = internalFaces_[i];
        zoneFaceU[facei] = Velocity;
    }

    // Included patches
    forAll(includedFaces_, patchi)
    {
        forAll(includedFaces_[patchi], i)
        {
            label patchFacei = includedFaces_[patchi][i];

            zoneFaceU.boundaryFieldRef()[patchi][patchFacei] = Velocity;
        }
    }

    // Excluded patches
    forAll(excludedFaces_, patchi)
    {
        forAll(excludedFaces_[patchi], i)
        {
            label patchFacei = excludedFaces_[patchi][i];

            zoneFaceU.boundaryFieldRef()[patchi][patchFacei] = Velocity;
        }
    }
}


void Foam::MRFTranslatingZone::correctBoundaryVelocity(volVectorField& U) const
{
    if (!active_)
    {
        return;
    }

    const vector Velocity = this->frameVelocity();

    // Included patches
  //  volVectorField::Boundary& Ubf = U.boundaryFieldRef();

    forAll(includedFaces_, patchi)
    {
        vectorField pfld(U.boundaryFieldRef()[patchi]);

        forAll(includedFaces_[patchi], i)
        {
            label patchFacei = includedFaces_[patchi][i];

            pfld[patchFacei] = Velocity;
        }

        U.boundaryFieldRef()[patchi] == pfld;
    }
}


void Foam::MRFTranslatingZone::writeData(Ostream& os) const
{
    os  << nl;
    os.beginBlock(name_);

    os.writeEntry("active", active_);
    os.writeEntry("cellZone", cellZoneName_);
    frameVel_->writeData(os);

    if (excludedPatchNames_.size())
    {
        os.writeEntry("nonMovingPatches", excludedPatchNames_);
    }

    os.endBlock();
}


bool Foam::MRFTranslatingZone::read(const dictionary& dict)
{
    coeffs_ = dict;

    active_ = coeffs_.getOrDefault("active", true);
    coeffs_.readEntry("cellZone", cellZoneName_);
    cellZoneID_ = mesh_.cellZones().findZoneID(cellZoneName_);

    return true;
}


void Foam::MRFTranslatingZone::update()
{
    if (mesh_.topoChanging())
    {
        setMRFFaces();
    }
}


// ************************************************************************* //
