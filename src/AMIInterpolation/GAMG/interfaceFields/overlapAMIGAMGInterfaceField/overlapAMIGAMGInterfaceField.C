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

#include "overlapAMIGAMGInterfaceField.H"

#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(overlapAMIGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        overlapAMIGAMGInterfaceField,
        lduInterface
    );
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        overlapAMIGAMGInterfaceField,
        lduInterfaceField
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::overlapAMIGAMGInterfaceField::overlapAMIGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    GAMGInterfaceField(GAMGCp, fineInterface),
    overlapAMIInterface_(refCast<const overlapAMIGAMGInterface>(GAMGCp)),
    doTransform_(false),
    rank_(0)
{
    const overlapAMILduInterfaceField& p =
        refCast<const overlapAMILduInterfaceField>(fineInterface);

    doTransform_ = p.doTransform();
    rank_ = p.rank();
}


Foam::overlapAMIGAMGInterfaceField::overlapAMIGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const bool doTransform,
    const int rank
)
:
    GAMGInterfaceField(GAMGCp, doTransform, rank),
    overlapAMIInterface_(refCast<const overlapAMIGAMGInterface>(GAMGCp)),
    doTransform_(doTransform),
    rank_(rank)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::overlapAMIGAMGInterfaceField::~overlapAMIGAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::overlapAMIGAMGInterfaceField::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
	const lduAddressing& lduAddr,
	const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    // Get neighbouring field

    const labelList& nbrFaceCells =
        lduAddr.patchAddr
        (
            overlapAMIInterface_.neighbPatchID()
        );

    solveScalarField pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    if (overlapAMIInterface_.owner())
    {
        pnf = overlapAMIInterface_.AMI().interpolateToSource(pnf);
    }
    else
    {
        pnf = overlapAMIInterface_.neighbPatch().AMI().interpolateToTarget(pnf);
    }

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}


// ************************************************************************* //
