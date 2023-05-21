/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvMatrix.H"
#include "cyclicACMIFvPatchField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicACMIFvPatchField<Type>::cyclicACMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicACMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF),
    cyclicACMIPatch_(refCast<const cyclicACMIFvPatch>(p))
{}


template<class Type>
Foam::cyclicACMIFvPatchField<Type>::cyclicACMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicACMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF, dict, dict.found("value")),
    cyclicACMIPatch_(refCast<const cyclicACMIFvPatch>(p, dict))
{
    if (!isA<cyclicACMIFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    if (!dict.found("value") && this->coupled())
    {
        // Extra check: make sure that the non-overlap patch is before
        // this so it has actually been read - evaluate will crash otherwise

        const GeometricField<Type, fvPatchField, volMesh>& fld =
            static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
            (
                this->primitiveField()
            );
        if (!fld.boundaryField().set(cyclicACMIPatch_.nonOverlapPatchID()))
        {
            FatalIOErrorInFunction(dict)
                << "    patch " << p.name()
                << " of field " << this->internalField().name()
                << " refers to non-overlap patch "
                << cyclicACMIPatch_.cyclicACMIPatch().nonOverlapPatchName()
                << " which is not constructed yet." << nl
                << "    Either supply an initial value or change the ordering"
                << " in the file"
                << exit(FatalIOError);
        }

        this->evaluate(Pstream::commsTypes::blocking);
    }
}


template<class Type>
Foam::cyclicACMIFvPatchField<Type>::cyclicACMIFvPatchField
(
    const cyclicACMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicACMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicACMIPatch_(refCast<const cyclicACMIFvPatch>(p))
{
    if (!isA<cyclicACMIFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }
}



template<class Type>
Foam::cyclicACMIFvPatchField<Type>::cyclicACMIFvPatchField
(
    const cyclicACMIFvPatchField<Type>& ptf
)
:
    cyclicACMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    cyclicACMIPatch_(ptf.cyclicACMIPatch_)
{}


template<class Type>
Foam::cyclicACMIFvPatchField<Type>::cyclicACMIFvPatchField
(
    const cyclicACMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicACMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    cyclicACMIPatch_(ptf.cyclicACMIPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::cyclicACMIFvPatchField<Type>::coupled() const
{
    return cyclicACMIPatch_.coupled();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicACMIFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->primitiveField();
    //const cyclicACMIPolyPatch& cpp = cyclicACMIPatch_.cyclicACMIPatch();

    // By pass polyPatch to get nbrId. Instead use cyclicAMIFvPatch virtual
    // neighbPatch()
    const cyclicACMIFvPatch& neighbPatch = cyclicACMIPatch_.neighbPatch();
    const labelUList& nbrFaceCells = neighbPatch.faceCells();

    tmp<Field<Type>> tpnf
    (
        cyclicACMIPatch_.interpolate
        (
            Field<Type>
            (
                iField,
                nbrFaceCells
                //cpp.neighbPatch().faceCells()
            )
        )
    );

    if (doTransform())
    {
        tpnf.ref() = transform(forwardT(), tpnf());
    }

    return tpnf;
}


template<class Type>
const Foam::cyclicACMIFvPatchField<Type>&
Foam::cyclicACMIFvPatchField<Type>::neighbourPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->primitiveField()
        );

    return refCast<const cyclicACMIFvPatchField<Type>>
    (
        fld.boundaryField()[cyclicACMIPatch_.neighbPatchID()]
    );
}


template<class Type>
const Foam::fvPatchField<Type>&
Foam::cyclicACMIFvPatchField<Type>::nonOverlapPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->primitiveField()
        );

    // WIP: Needs to re-direct nonOverlapPatchID to new patchId for assembly?
    return fld.boundaryField()[cyclicACMIPatch_.nonOverlapPatchID()];
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::updateInterfaceMatrix
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
    // note: only applying coupled contribution

//     const labelUList& nbrFaceCellsCoupled =
//         lduAddr.patchAddr
//         (
//             cyclicACMIPatch_.cyclicACMIPatch().neighbPatchID()
//         );

    const labelUList& nbrFaceCellsCoupled =
        lduAddr.patchAddr(cyclicACMIPatch_.neighbPatchID());

    solveScalarField pnf(psiInternal, nbrFaceCellsCoupled);

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    pnf = cyclicACMIPatch_.interpolate(pnf);

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{
    // note: only applying coupled contribution

    const labelUList& nbrFaceCellsCoupled =
        lduAddr.patchAddr(cyclicACMIPatch_.neighbPatchID());

    Field<Type> pnf(psiInternal, nbrFaceCellsCoupled);

    // Transform according to the transformation tensors
    transformCoupleField(pnf);

    pnf = cyclicACMIPatch_.interpolate(pnf);

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix
)
{
    const scalarField& mask = cyclicACMIPatch_.cyclicACMIPatch().mask();

    // Nothing to be done by the AMI, but re-direct to non-overlap patch
    // with non-overlap patch weights
    const fvPatchField<Type>& npf = nonOverlapPatchField();

    const_cast<fvPatchField<Type>&>(npf).manipulateMatrix(matrix, 1.0 - mask);
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix,
    const label mat,
    const direction cmpt
)
{
    if (this->cyclicACMIPatch().owner())
    {
        label index = this->patch().index();

        const label globalPatchID =
            matrix.lduMeshAssembly().patchLocalToGlobalMap()[mat][index];

        const Field<scalar> intCoeffsCmpt
        (
            matrix.internalCoeffs()[globalPatchID].component(cmpt)
        );

        const Field<scalar> boundCoeffsCmpt
        (
            matrix.boundaryCoeffs()[globalPatchID].component(cmpt)
        );

        tmp<Field<scalar>> tintCoeffs(coeffs(matrix, intCoeffsCmpt, mat));
        tmp<Field<scalar>> tbndCoeffs(coeffs(matrix, boundCoeffsCmpt, mat));
        const Field<scalar>& intCoeffs = tintCoeffs.ref();
        const Field<scalar>& bndCoeffs = tbndCoeffs.ref();

        const labelUList& u = matrix.lduAddr().upperAddr();
        const labelUList& l = matrix.lduAddr().lowerAddr();

        label subFaceI = 0;

        const labelList& faceMap =
            matrix.lduMeshAssembly().faceBoundMap()[mat][index];

        forAll (faceMap, j)
        {
            label globalFaceI = faceMap[j];

            const scalar boundCorr = -bndCoeffs[subFaceI];
            const scalar intCorr = -intCoeffs[subFaceI];

            matrix.upper()[globalFaceI] += boundCorr;
            matrix.diag()[u[globalFaceI]] -= intCorr;
            matrix.diag()[l[globalFaceI]] -= boundCorr;

            if (matrix.asymmetric())
            {
                matrix.lower()[globalFaceI] += intCorr;
            }
            subFaceI++;
        }

        // Set internalCoeffs and boundaryCoeffs in the assembly matrix
        // on clyclicAMI patches to be used in the individual matrix by
        // matrix.flux()
        if (matrix.psi(mat).mesh().fluxRequired(this->internalField().name()))
        {
            matrix.internalCoeffs().set
            (
                globalPatchID, intCoeffs*pTraits<Type>::one
            );
            matrix.boundaryCoeffs().set
            (
                globalPatchID, bndCoeffs*pTraits<Type>::one
            );

            const label nbrPathID =
                cyclicACMIPatch_.cyclicACMIPatch().neighbPatchID();

            const label nbrGlobalPatchID =
                matrix.lduMeshAssembly().patchLocalToGlobalMap()
                [mat][nbrPathID];

            matrix.internalCoeffs().set
            (
                nbrGlobalPatchID, intCoeffs*pTraits<Type>::one
            );
            matrix.boundaryCoeffs().set
            (
                nbrGlobalPatchID, bndCoeffs*pTraits<Type>::one
            );
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Foam::scalar>>
Foam::cyclicACMIFvPatchField<Type>::coeffs
(
    fvMatrix<Type>& matrix,
    const Field<scalar>& coeffs,
    const label mat
) const
{
    const label index(this->patch().index());

    const label nSubFaces
    (
        matrix.lduMeshAssembly().cellBoundMap()[mat][index].size()
    );

    Field<scalar> mapCoeffs(nSubFaces, Zero);

    const scalarListList& srcWeight =
        cyclicACMIPatch_.cyclicACMIPatch().AMI().srcWeights();

    const scalarField& mask = cyclicACMIPatch_.cyclicACMIPatch().mask();

    const scalar tol = cyclicACMIPolyPatch::tolerance();
    label subFaceI = 0;
    forAll(mask, faceI)
    {
        const scalarList& w = srcWeight[faceI];
        for(label i=0; i<w.size(); i++)
        {
            if (mask[faceI] > tol)
            {
                const label localFaceId =
                    matrix.lduMeshAssembly().facePatchFaceMap()
                    [mat][index][subFaceI];
                mapCoeffs[subFaceI] = w[i]*coeffs[localFaceId];
            }
            subFaceI++;
        }
    }

    return tmp<Field<scalar>>(new Field<scalar>(mapCoeffs));
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::updateCoeffs()
{
    // Update non-overlap patch - some will implement updateCoeffs, and
    // others will implement evaluate

    // Pass in (1 - mask) to give non-overlap patch the chance to do
    // manipulation of non-face based data

    const scalarField& mask = cyclicACMIPatch_.cyclicACMIPatch().mask();
    const fvPatchField<Type>& npf = nonOverlapPatchField();
    const_cast<fvPatchField<Type>&>(npf).updateWeightedCoeffs(1.0 - mask);
}


template<class Type>
void Foam::cyclicACMIFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
