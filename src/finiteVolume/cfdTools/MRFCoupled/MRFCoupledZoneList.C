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

#include "MRFCoupledZoneList.H"

#include "volFields.H"
#include "fixedValueFvsPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFCoupledZoneList::MRFCoupledZoneList
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    PtrList<MRFCoupledZone>(),
    mesh_(mesh)
{
    reset(dict);

    active(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::MRFCoupledZoneList::~MRFCoupledZoneList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::MRFCoupledZoneList::omega() const
{
    tmp<volVectorField> tMRFZonesOmega
    (
        new volVectorField
        (
            IOobject
            (
                "MRFZonesOmega",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimless/dimTime, vector::zero)
        )
    );
    volVectorField& MRFZonesOmega = tMRFZonesOmega.ref();

    forAll (*this, i)
    {
        operator[](i).addOmega(MRFZonesOmega);
    }

    return tMRFZonesOmega;
}


Foam::tmp<Foam::volVectorField> Foam::MRFCoupledZoneList::origin() const
{
    tmp<volVectorField> tOrigin
    (
        new volVectorField
        (
            IOobject
            (
                "MRFZonesOrigin",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimless/dimTime, vector::zero)
        )
    );
    volVectorField& MRFZonesOrigin = tOrigin.ref();

    forAll (*this, i)
    {
        operator[](i).addOrigin(MRFZonesOrigin);
    }

    return tOrigin;
}


Foam::tmp<Foam::surfaceVectorField> Foam::MRFCoupledZoneList::faceU() const
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


bool Foam::MRFCoupledZoneList::active(const bool warn) const
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


void Foam::MRFCoupledZoneList::reset(const dictionary& dict)
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
                new MRFCoupledZone(name, mesh_, modelDict)
            );
        }
    }
}


bool Foam::MRFCoupledZoneList::read(const dictionary& dict)
{
    bool allOk = true;
    forAll(*this, i)
    {
        MRFCoupledZone& pm = this->operator[](i);
        bool ok = pm.read(dict.subDict(pm.name()));
        allOk = (allOk && ok);
    }
    return allOk;
}


bool Foam::MRFCoupledZoneList::writeData(Ostream& os) const
{
    forAll(*this, i)
    {
        os  << nl;
        this->operator[](i).writeData(os);
    }

    return os.good();
}


void Foam::MRFCoupledZoneList::addAcceleration
(
    const volVectorField& U,
    volVectorField& ddtU
) const
{
    forAll(*this, i)
     {
        operator[](i).addCoriolis(U, ddtU);
    }
}


void Foam::MRFCoupledZoneList::addAcceleration(fvVectorMatrix& UEqn) const
{
    forAll(*this, i)
    {
        operator[](i).addCoriolis(UEqn);
    }
}


void Foam::MRFCoupledZoneList::addAcceleration
(
    const volScalarField& rho,
    fvVectorMatrix& UEqn
) const
{
    forAll(*this, i)
    {
        operator[](i).addCoriolis(rho, UEqn);
    }
}


Foam::tmp<Foam::volVectorField> Foam::MRFCoupledZoneList::DDt
(
    const volVectorField& U
) const
{
    tmp<volVectorField> tacceleration
    (
        new volVectorField
        (
            IOobject
            (
                "MRFCoupledZoneList:acceleration",
                U.mesh().time().timeName(),
                U.mesh()
            ),
            U.mesh(),
            dimensionedVector(U.dimensions()/dimTime, Zero)
        )
    );
    volVectorField& acceleration = tacceleration.ref();

    forAll(*this, i)
    {
        operator[](i).addCoriolis(U, acceleration);
    }

    return tacceleration;
}


Foam::tmp<Foam::volVectorField> Foam::MRFCoupledZoneList::DDt
(
    const volScalarField& rho,
    const volVectorField& U
) const
{
    return rho*DDt(U);
}


void Foam::MRFCoupledZoneList::makeRelative(volVectorField& U) const
{
    forAll(*this, i)
    {
        operator[](i).makeRelative(U);
    }
}


void Foam::MRFCoupledZoneList::makeRelative(surfaceScalarField& phi) const
{
    forAll(*this, i)
    {
        operator[](i).makeRelative(phi);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::MRFCoupledZoneList::relative
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
Foam::MRFCoupledZoneList::relative
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
Foam::MRFCoupledZoneList::relative
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


void Foam::MRFCoupledZoneList::makeRelative
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


void Foam::MRFCoupledZoneList::makeAbsolute(volVectorField& U) const
{
    forAll(*this, i)
    {
        operator[](i).makeAbsolute(U);
    }
}


void Foam::MRFCoupledZoneList::makeAbsolute(surfaceScalarField& phi) const
{
    forAll(*this, i)
    {
        operator[](i).makeAbsolute(phi);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::MRFCoupledZoneList::absolute
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


void Foam::MRFCoupledZoneList::makeAbsolute
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


void Foam::MRFCoupledZoneList::correctBoundaryVelocity(volVectorField& U) const
{
    forAll(*this, i)
    {
        operator[](i).correctBoundaryVelocity(U);
    }
}


void Foam::MRFCoupledZoneList::correctBoundaryFlux
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


void Foam::MRFCoupledZoneList::update()
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
    const MRFCoupledZoneList& models
)
{
    models.writeData(os);
    return os;
}


// ************************************************************************* //
