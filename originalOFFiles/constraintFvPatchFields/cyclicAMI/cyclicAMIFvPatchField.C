/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
#include "volFields.H"
//#include "cylicFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p))
{}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(p, iF, dict, dict.found("value")),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p, dict))
{
    if (!isA<cyclicAMIFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    if (!dict.found("value"))
    {
        if (this->coupled())
        {
            this->evaluate(Pstream::commsTypes::blocking);
        }
        else
        {
            fvPatchField<Type>::operator=(this->patchInternalField());
        }
    }
}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p))
{
    if (!isA<cyclicAMIFvPatch>(this->patch()))
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
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_)
{}


template<class Type>
Foam::cyclicAMIFvPatchField<Type>::cyclicAMIFvPatchField
(
    const cyclicAMIFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMILduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::cyclicAMIFvPatchField<Type>::coupled() const
{
    return cyclicAMIPatch_.coupled();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->primitiveField();

    // By pass polyPatch to get nbrId. Instead use cyclicAMIFvPatch virtual
    // neighbPatch()
    const cyclicAMIFvPatch& neighbPatch = cyclicAMIPatch_.neighbPatch();
    const labelUList& nbrFaceCells = neighbPatch.faceCells();

    Field<Type> pnf(iField, nbrFaceCells);

    const word& fieldName = this->internalField().name();

    tmp<Field<Type>> tpnf;
    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        Field<Type> pnfInternal(iField, cyclicAMIPatch_.faceCells());

        tpnf = cyclicAMIPatch_.interpolate(pnf, pnfInternal);
    }
    else
    {
        tpnf = cyclicAMIPatch_.interpolate(pnf);
    }

    if (doTransform())
    {
		if (fieldName.startsWith("U"))
		{
			if (fieldName == "U")
			{
				tpnf.ref() = transform(forwardT(), tpnf());
			}
			else
			{
				const fvMesh& mesh = this->patch().boundaryMesh().mesh();
				const vectorField& Ui = mesh.objectRegistry::lookupObject<volVectorField>("U").primitiveField();
				vectorField Up(Ui, nbrFaceCells);
				vectorField UpInterp = cyclicAMIPatch_.interpolate(Up);

				UpInterp = transform(forwardT(), UpInterp);

				if (fieldName == "U.component(0)")
				{
					tpnf.ref() = UpInterp.component(0)*pTraits<Type>::one;
				}
				else if (fieldName == "U.component(1)")
				{
					tpnf.ref() = UpInterp.component(1)*pTraits<Type>::one;
				}
				else if  (fieldName == "U.component(2)")
				{
					tpnf.ref() = UpInterp.component(2)*pTraits<Type>::one;
				}
			}
		}
		else
		{
			tpnf.ref() = transform(forwardT(), tpnf());
		}
    }

    return tpnf;
}


template<class Type>
const Foam::cyclicAMIFvPatchField<Type>&
Foam::cyclicAMIFvPatchField<Type>::neighbourPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
        static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->primitiveField()
        );

    return refCast<const cyclicAMIFvPatchField<Type>>
    (
        fld.boundaryField()[cyclicAMIPatch_.neighbPatchID()]
    );
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
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
    const labelUList& nbrFaceCells =
        lduAddr.patchAddr(cyclicAMIPatch_.neighbPatchID());

    solveScalarField pnf(psiInternal, nbrFaceCells);

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        solveScalarField pif(psiInternal, faceCells);
        pnf = cyclicAMIPatch_.interpolate(pnf, pif);
    }
    else
    {
        pnf = cyclicAMIPatch_.interpolate(pnf);
    }

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
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
    const labelUList& nbrFaceCells =
        lduAddr.patchAddr(cyclicAMIPatch_.neighbPatchID());

    Field<Type> pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf);

    if (cyclicAMIPatch_.applyLowWeightCorrection())
    {
        Field<Type> pif(psiInternal, cyclicAMIPatch_.faceCells());
        pnf = cyclicAMIPatch_.interpolate(pnf, pif);
    }
    else
    {
        pnf = cyclicAMIPatch_.interpolate(pnf);
    }

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, faceCells, coeffs, pnf);
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix,
    const label mat,
    const direction cmpt
)
{

    if (this->cyclicAMIPatch().owner())
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
                cyclicAMIPatch_.cyclicAMIPatch().neighbPatchID();

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
Foam::cyclicAMIFvPatchField<Type>::coeffs
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
        cyclicAMIPatch_.cyclicAMIPatch().AMI().srcWeights();

    label subFaceI = 0;
    forAll(*this, faceI)
    {
        const scalarList& w = srcWeight[faceI];
        for(label i=0; i<w.size(); i++)
        {
            const label localFaceId =
                matrix.lduMeshAssembly().facePatchFaceMap()[mat][index][subFaceI];
            mapCoeffs[subFaceI] = w[i]*coeffs[localFaceId];
            subFaceI++;
        }
    }

    return tmp<Field<scalar>>(new Field<scalar>(mapCoeffs));
}


template<class Type>
template<class Type2>
void Foam::cyclicAMIFvPatchField<Type>::collectStencilData
(
    const refPtr<mapDistribute>& mapPtr,
    const labelListList& stencil,
    const Type2& data,
    List<Type2>& expandedData
)
{
    expandedData.setSize(stencil.size());
    if (mapPtr.valid())
    {
        Type2 work(data);
        mapPtr().distribute(work);

        forAll(stencil, facei)
        {
            const labelList& slots = stencil[facei];
            expandedData[facei].append
            (
                UIndirectList<typename Type2::value_type>(work, slots)
            );
        }
    }
    else
    {
        forAll(stencil, facei)
        {
            const labelList& slots = stencil[facei];
            expandedData[facei].append
            (
                UIndirectList<typename Type2::value_type>(data, slots)
            );
        }
    }
}


template<class Type>
void Foam::cyclicAMIFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //

