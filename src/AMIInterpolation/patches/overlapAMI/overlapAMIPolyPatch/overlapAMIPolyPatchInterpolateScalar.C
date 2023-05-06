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

#include "overlapAMIPolyPatch.H"
#include "psiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::scalarField >
Foam::overlapAMIPolyPatch::expandData(const scalarField& pf, label cmpt) const
{
    // Check and expand the field from patch size to zone size
    if (pf.size() != this->size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > overlapAMIPolyPatch::expandData"
            "("
            "    const Field<Type>& pf"
            ") const"
        )   << "Incorrect patch field size.  Field size: "
            << pf.size() << " patch size: " << this->size()
            << abort(FatalError);
    }

    const label ncp = nCopies();

	const scalar myAngle = 360.0/scalar(ncp);

    tmp<scalarField > texpandField
    (
	   new scalarField(ncp*pf.size(), 0.0)
    );

    scalarField& expandField = texpandField.ref();

    const labelUList& faceCells = this->faceCells();
    const vectorField& Ui = this->boundaryMesh().mesh().objectRegistry::lookupObject<volVectorField>("U").primitiveField();
    vectorField Up(Ui, faceCells);

    for (label copyI = 0; copyI < ncp; copyI++)
    {
    	// Calculate transform
		const tensor curRotation = this->RodriguesRotation(rotationAxis_, copyI*myAngle);

		const label offset = copyI*pf.size();

		forAll (pf, faceI)
		{
			 const label zId = this->whichFace(this->start() + faceI);
			 vector UTransf = Foam::transform(curRotation, Up[faceI]);
			 expandField[offset + zId] = UTransf.component(cmpt);
		}
    }



    return texpandField;
}

// ************************************************************************* //
