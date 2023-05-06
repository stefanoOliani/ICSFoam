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

#include "phaseLagAMIFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "transform.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseLagAMIFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, phaseLagAMIFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::phaseLagAMIFvPatch::coupled() const
{
    return Pstream::parRun() || (this->size() && neighbFvPatch().size());
}


void Foam::phaseLagAMIFvPatch::makeWeights(scalarField& w) const
{
    if (coupled())
    {
        const phaseLagAMIFvPatch& nbrPatch = neighbFvPatch();

        const scalarField deltas(nf() & coupledFvPatch::delta());

        tmp<scalarField> tnbrDeltas;
        if (applyLowWeightCorrection())
        {
            tnbrDeltas =
                interpolate
                (
                    nbrPatch.nf() & nbrPatch.coupledFvPatch::delta(),
                    scalarField(this->size(), 1.0)
                );
        }
        else
        {
            tnbrDeltas =
                interpolate(nbrPatch.nf() & nbrPatch.coupledFvPatch::delta());
        }

        const scalarField& nbrDeltas = tnbrDeltas();

        forAll(deltas, facei)
        {
            // Note use of mag
            scalar di = mag(deltas[facei]);
            scalar dni = mag(nbrDeltas[facei]);

            w[facei] = dni/(di + dni);
        }
    }
    else
    {
        // Behave as uncoupled patch
        fvPatch::makeWeights(w);
    }
}


void Foam::phaseLagAMIFvPatch::makeDeltaCoeffs(scalarField& coeffs) const
{
    // Apply correction to default coeffs
}


void Foam::phaseLagAMIFvPatch::makeNonOrthoDeltaCoeffs(scalarField& coeffs) const
{
    // Apply correction to default coeffs
    //coeffs = Zero;
}


void Foam::phaseLagAMIFvPatch::makeNonOrthoCorrVectors(vectorField& vecs) const
{
    // Apply correction to default vectors
    //vecs = Zero;
}


Foam::tmp<Foam::vectorField> Foam::phaseLagAMIFvPatch::delta() const
{
    const phaseLagAMIFvPatch& nbrPatch = neighbFvPatch();

    if (coupled())
    {
        const vectorField patchD(coupledFvPatch::delta());

        tmp<vectorField> tnbrPatchD;
        if (applyLowWeightCorrection())
        {
            tnbrPatchD =
                interpolate
                (
                    nbrPatch.coupledFvPatch::delta(),
                    vectorField(this->size(), Zero)
                );
        }
        else
        {
            tnbrPatchD = interpolate(nbrPatch.coupledFvPatch::delta());
        }

        const vectorField& nbrPatchD = tnbrPatchD();

        auto tpdv = tmp<vectorField>::New(patchD.size());
        vectorField& pdv = tpdv.ref();

        // do the transformation if necessary
        if (parallel())
        {
            forAll(patchD, facei)
            {
                const vector& ddi = patchD[facei];
                const vector& dni = nbrPatchD[facei];

                pdv[facei] = ddi - dni;
            }
        }
        else
        {
            forAll(patchD, facei)
            {
                const vector& ddi = patchD[facei];
                const vector& dni = nbrPatchD[facei];

                pdv[facei] = ddi - transform(forwardT()[0], dni);
            }
        }

        return tpdv;
    }
    else
    {
        return coupledFvPatch::delta();
    }
}


Foam::tmp<Foam::labelField> Foam::phaseLagAMIFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::phaseLagAMIFvPatch::interfaceInternalField
(
    const labelUList& internalData,
	const labelUList& faceCells
) const
{
    return patchInternalField(internalData, faceCells);
}


Foam::tmp<Foam::labelField> Foam::phaseLagAMIFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    return neighbFvPatch().patchInternalField(iF);
}


void Foam::phaseLagAMIFvPatch::movePoints()
{
}


// ************************************************************************* //
