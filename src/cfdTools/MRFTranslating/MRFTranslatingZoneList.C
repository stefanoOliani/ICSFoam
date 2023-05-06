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

#include "MRFTranslatingZoneList.H"

#include "volFields.H"
#include "fixedValueFvsPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFTranslatingZoneList::MRFTranslatingZoneList
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    PtrList<MRFTranslatingZone>(),
    mesh_(mesh)
{
    reset(dict);

    active(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::MRFTranslatingZoneList::~MRFTranslatingZoneList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::MRFTranslatingZoneList::frameVel() const
{
    tmp<volVectorField> tMRFZonesFrameVel
    (
        new volVectorField
        (
            IOobject
            (
                "MRFZonesFrameVel",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimless/dimTime, vector::zero)
        )
    );
    volVectorField& MRFZonesFrameVel = tMRFZonesFrameVel.ref();

    forAll (*this, i)
    {
        operator[](i).addFrameVelocity(MRFZonesFrameVel);
    }

    return tMRFZonesFrameVel;
}



Foam::tmp<Foam::surfaceVectorField> Foam::MRFTranslatingZoneList::faceU() const
{
    tmp<surfaceVectorField> tMRFZonesFaceU
    (
        new surfaceVectorField
        (
            IOobject
            (
                "MRFZonesFaceU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimVelocity, vector::zero)
        )
    );
    surfaceVectorField& MRFZonesFaceU = tMRFZonesFaceU.ref();

    forAll(*this, i)
    {
        operator[](i).faceU(MRFZonesFaceU);
    }

    return tMRFZonesFaceU;
}


bool Foam::MRFTranslatingZoneList::active(const bool warn) const
{
    bool a = false;
    forAll(*this, i)
    {
        a = a || this->operator[](i).active();
    }

    if (warn && this->size() && !a)
    {
        Info<< "    No MRF zones active" << endl;
    }

    return a;
}


void Foam::MRFTranslatingZoneList::reset(const dictionary& dict)
{
    label count = 0;
    for (const entry& dEntry : dict)
    {
        if (dEntry.isDict())
        {
            ++count;
        }
    }

    this->resize(count);

    count = 0;
    for (const entry& dEntry : dict)
    {
        if (dEntry.isDict())
        {
            const word& name = dEntry.keyword();
            const dictionary& modelDict = dEntry.dict();

            Info<< "    creating MRF zone: " << name << endl;

            this->set
            (
                count++,
                new MRFTranslatingZone(name, mesh_, modelDict)
            );
        }
    }
}


bool Foam::MRFTranslatingZoneList::read(const dictionary& dict)
{
    bool allOk = true;
    forAll(*this, i)
    {
        MRFTranslatingZone& pm = this->operator[](i);
        bool ok = pm.read(dict.subDict(pm.name()));
        allOk = (allOk && ok);
    }
    return allOk;
}


bool Foam::MRFTranslatingZoneList::writeData(Ostream& os) const
{
    forAll(*this, i)
    {
        os  << nl;
        this->operator[](i).writeData(os);
    }

    return os.good();
}


void Foam::MRFTranslatingZoneList::makeRelative(volVectorField& U) const
{
    forAll(*this, i)
    {
        operator[](i).makeRelative(U);
    }
}


void Foam::MRFTranslatingZoneList::makeRelative(surfaceScalarField& phi) const
{
    forAll(*this, i)
    {
        operator[](i).makeRelative(phi);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::MRFTranslatingZoneList::relative
(
    const tmp<surfaceScalarField>& tphi
) const
{
    if (size())
    {
        tmp<surfaceScalarField> rphi
        (
            New
            (
                tphi,
                "relative(" + tphi().name() + ')',
                tphi().dimensions(),
                true
            )
        );

        makeRelative(rphi.ref());

        tphi.clear();

        return rphi;
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


Foam::tmp<Foam::FieldField<Foam::fvsPatchField, Foam::scalar>>
Foam::MRFTranslatingZoneList::relative
(
    const tmp<FieldField<fvsPatchField, scalar>>& tphi
) const
{
    if (size())
    {
        tmp<FieldField<fvsPatchField, scalar>> rphi(New(tphi, true));

        forAll(*this, i)
        {
            operator[](i).makeRelative(rphi.ref());
        }

        tphi.clear();

        return rphi;
    }
    else
    {
        return tmp<FieldField<fvsPatchField, scalar>>(tphi, true);
    }
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::MRFTranslatingZoneList::relative
(
    const tmp<Field<scalar>>& tphi,
    const label patchi
) const
{
    if (size())
    {
        tmp<Field<scalar>> rphi(New(tphi, true));

        forAll(*this, i)
        {
            operator[](i).makeRelative(rphi.ref(), patchi);
        }

        tphi.clear();

        return rphi;
    }
    else
    {
        return tmp<Field<scalar>>(tphi, true);
    }
}


void Foam::MRFTranslatingZoneList::makeRelative
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    forAll(*this, i)
    {
        operator[](i).makeRelative(rho, phi);
    }
}


void Foam::MRFTranslatingZoneList::makeAbsolute(volVectorField& U) const
{
    forAll(*this, i)
    {
        operator[](i).makeAbsolute(U);
    }
}


void Foam::MRFTranslatingZoneList::makeAbsolute(surfaceScalarField& phi) const
{
    forAll(*this, i)
    {
        operator[](i).makeAbsolute(phi);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::MRFTranslatingZoneList::absolute
(
    const tmp<surfaceScalarField>& tphi
) const
{
    if (size())
    {
        tmp<surfaceScalarField> rphi
        (
            New
            (
                tphi,
                "absolute(" + tphi().name() + ')',
                tphi().dimensions(),
                true
            )
        );

        makeAbsolute(rphi.ref());

        tphi.clear();

        return rphi;
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


void Foam::MRFTranslatingZoneList::makeAbsolute
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    forAll(*this, i)
    {
        operator[](i).makeAbsolute(rho, phi);
    }
}


void Foam::MRFTranslatingZoneList::correctBoundaryVelocity(volVectorField& U) const
{
    forAll(*this, i)
    {
        operator[](i).correctBoundaryVelocity(U);
    }
}


void Foam::MRFTranslatingZoneList::correctBoundaryFlux
(
    const volVectorField& U,
    surfaceScalarField& phi
) const
{
    FieldField<fvsPatchField, scalar> Uf
    (
        relative(mesh_.Sf().boundaryField() & U.boundaryField())
    );


    surfaceScalarField::Boundary& phibf = phi.boundaryFieldRef();

    forAll(mesh_.boundary(), patchi)
    {
        if
        (
            isA<fixedValueFvsPatchScalarField>(phibf[patchi])
        )
        {
            phibf[patchi] == Uf[patchi];
        }
    }
}


void Foam::MRFTranslatingZoneList::update()
{
    if (mesh_.topoChanging())
    {
        forAll(*this, i)
        {
            operator[](i).update();
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const MRFTranslatingZoneList& models
)
{
    models.writeData(os);
    return os;
}


// ************************************************************************* //
