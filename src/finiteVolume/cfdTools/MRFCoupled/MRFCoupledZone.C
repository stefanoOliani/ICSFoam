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

#include "MRFCoupledZone.H"

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
    defineTypeNameAndDebug(MRFCoupledZone, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::MRFCoupledZone::setMRFFaces()
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

Foam::MRFCoupledZone::MRFCoupledZone
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
        coeffs_.getOrDefault<wordRes>("nonRotatingPatches", wordRes())
    ),
    origin_(coeffs_.get<vector>("origin")),
    axis_(coeffs_.get<vector>("axis").normalise()),
    omega_(Function1<scalar>::New("omega", coeffs_))
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

Foam::vector Foam::MRFCoupledZone::Omega() const
{
    return omega_->value(mesh_.time().timeOutputValue())*axis_;
}


void Foam::MRFCoupledZone::addOmega(volVectorField& omg) const
{
    const vector rotVel = Omega();

    // Set omega in all cells of the rotating cell zone
    const labelList& cells = mesh_.cellZones()[cellZoneID_];

    forAll (cells, i)
    {
        omg[cells[i]] = rotVel;
    }
}

void Foam::MRFCoupledZone::addOrigin(volVectorField& orig) const
{
    // Set omega in all cells of the rotating cell zone
    const labelList& cells = mesh_.cellZones()[cellZoneID_];

    forAll (cells, i)
    {
        orig[cells[i]] = origin_;
    }
}


void Foam::MRFCoupledZone::addCoriolis
(
    const volVectorField& U,
    volVectorField& ddtU
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    vectorField& ddtUc = ddtU.primitiveFieldRef();
    const vectorField& Uc = U;

    const vector Omega = this->Omega();

    forAll(cells, i)
    {
        label celli = cells[i];
        ddtUc[celli] += (Omega ^ Uc[celli]);
    }
}


void Foam::MRFCoupledZone::addCoriolis(fvVectorMatrix& UEqn, const bool rhs) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const scalarField& V = mesh_.V();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();

    const vector Omega = this->Omega();

    if (rhs)
    {
        forAll(cells, i)
        {
            label celli = cells[i];
            Usource[celli] += V[celli]*(Omega ^ U[celli]);
        }
    }
    else
    {
        forAll(cells, i)
        {
            label celli = cells[i];
            Usource[celli] -= V[celli]*(Omega ^ U[celli]);
        }
    }
}


void Foam::MRFCoupledZone::addCoriolis
(
    const volScalarField& rho,
    fvVectorMatrix& UEqn,
    const bool rhs
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const scalarField& V = mesh_.V();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();

    const vector Omega = this->Omega();

    if (rhs)
    {
        forAll(cells, i)
        {
            label celli = cells[i];
            Usource[celli] += V[celli]*rho[celli]*(Omega ^ U[celli]);
        }
    }
    else
    {
        forAll(cells, i)
        {
            label celli = cells[i];
            Usource[celli] -= V[celli]*rho[celli]*(Omega ^ U[celli]);
        }
    }
}


void Foam::MRFCoupledZone::makeRelative(volVectorField& U) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const volVectorField& C = mesh_.C();

    const vector Omega = this->Omega();

    const labelList& cells = mesh_.cellZones()[cellZoneID_];

    forAll(cells, i)
    {
        label celli = cells[i];
        U[celli] -= (Omega ^ (C[celli] - origin_));
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
            Ubf[patchi][patchFacei] -=
                (Omega
              ^ (C.boundaryField()[patchi][patchFacei] - origin_));
        }
    }
}


void Foam::MRFCoupledZone::makeRelative(surfaceScalarField& phi) const
{
    makeRelativeRhoFlux(geometricOneField(), phi);
}


void Foam::MRFCoupledZone::makeRelative(FieldField<fvsPatchField, scalar>& phi) const
{
    makeRelativeRhoFlux(oneFieldField(), phi);
}


void Foam::MRFCoupledZone::makeRelative(Field<scalar>& phi, const label patchi) const
{
    makeRelativeRhoFlux(oneField(), phi, patchi);
}


void Foam::MRFCoupledZone::makeRelative
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    makeRelativeRhoFlux(rho, phi);
}


void Foam::MRFCoupledZone::makeAbsolute(volVectorField& U) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const volVectorField& C = mesh_.C();

    const vector Omega = this->Omega();

    const labelList& cells = mesh_.cellZones()[cellZoneID_];

    forAll(cells, i)
    {
        label celli = cells[i];
        U[celli] += (Omega ^ (C[celli] - origin_));
    }

    // Included patches
    volVectorField::Boundary& Ubf = U.boundaryFieldRef();

    forAll(includedFaces_, patchi)
    {
        forAll(includedFaces_[patchi], i)
        {
            label patchFacei = includedFaces_[patchi][i];
            Ubf[patchi][patchFacei] =
                (Omega ^ (C.boundaryField()[patchi][patchFacei] - origin_));
        }
    }

    // Excluded patches
    forAll(excludedFaces_, patchi)
    {
        forAll(excludedFaces_[patchi], i)
        {
            label patchFacei = excludedFaces_[patchi][i];
            Ubf[patchi][patchFacei] +=
                (Omega ^ (C.boundaryField()[patchi][patchFacei] - origin_));
        }
    }
}


void Foam::MRFCoupledZone::makeAbsolute(surfaceScalarField& phi) const
{
    makeAbsoluteRhoFlux(geometricOneField(), phi);
}


void Foam::MRFCoupledZone::makeAbsolute
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    makeAbsoluteRhoFlux(rho, phi);
}


void Foam::MRFCoupledZone::faceU
(
    surfaceVectorField& zoneFaceU
) const
{
    const surfaceVectorField& Cf = mesh_.Cf();

    const vector& Omega = this->Omega();

    // Internal faces
    forAll(internalFaces_, i)
    {
        label facei = internalFaces_[i];
        zoneFaceU[facei] = (Omega ^ (Cf[facei] - origin_));
    }

    // Included patches
    forAll(includedFaces_, patchi)
    {
        forAll(includedFaces_[patchi], i)
        {
            label patchFacei = includedFaces_[patchi][i];

            zoneFaceU.boundaryFieldRef()[patchi][patchFacei] =
                (Omega ^ (Cf.boundaryField()[patchi][patchFacei] - origin_));
        }
    }

    // Excluded patches
    forAll(excludedFaces_, patchi)
    {
        forAll(excludedFaces_[patchi], i)
        {
            label patchFacei = excludedFaces_[patchi][i];

            zoneFaceU.boundaryFieldRef()[patchi][patchFacei] =
                (Omega ^ (Cf.boundaryField()[patchi][patchFacei] - origin_));
        }
    }
}


void Foam::MRFCoupledZone::correctBoundaryVelocity(volVectorField& U) const
{
    if (!active_)
    {
        return;
    }

    const vector Omega = this->Omega();

    // Included patches

    forAll(includedFaces_, patchi)
    {
        const vectorField& patchC = mesh_.Cf().boundaryField()[patchi];

        vectorField pfld(U.boundaryFieldRef()[patchi]);

        forAll(includedFaces_[patchi], i)
        {
            label patchFacei = includedFaces_[patchi][i];

            pfld[patchFacei] = (Omega ^ (patchC[patchFacei] - origin_));
        }

        U.boundaryFieldRef()[patchi] == pfld;
    }
}


void Foam::MRFCoupledZone::writeData(Ostream& os) const
{
    os  << nl;
    os.beginBlock(name_);

    os.writeEntry("active", active_);
    os.writeEntry("cellZone", cellZoneName_);
    os.writeEntry("origin", origin_);
    os.writeEntry("axis", axis_);
    omega_->writeData(os);

    if (excludedPatchNames_.size())
    {
        os.writeEntry("nonRotatingPatches", excludedPatchNames_);
    }

    os.endBlock();
}


bool Foam::MRFCoupledZone::read(const dictionary& dict)
{
    coeffs_ = dict;

    active_ = coeffs_.getOrDefault("active", true);
    coeffs_.readEntry("cellZone", cellZoneName_);
    cellZoneID_ = mesh_.cellZones().findZoneID(cellZoneName_);

    return true;
}


void Foam::MRFCoupledZone::update()
{
    if (mesh_.topoChanging())
    {
        setMRFFaces();
    }
}


// ************************************************************************* //
