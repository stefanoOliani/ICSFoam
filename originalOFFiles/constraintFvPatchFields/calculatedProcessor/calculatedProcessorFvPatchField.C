/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "calculatedProcessorFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::calculatedProcessorFvPatchField<Type>::calculatedProcessorFvPatchField
(
    const lduInterface& interface,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    procInterface_(refCast<const lduPrimitiveProcessorInterface>(interface)),
    sendBuf_(interface.faceCells().size()),
    receiveBuf_(interface.faceCells().size()),
    scalarSendBuf_(interface.faceCells().size()),
    scalarReceiveBuf_(interface.faceCells().size()),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{}


template<class Type>
Foam::calculatedProcessorFvPatchField<Type>::calculatedProcessorFvPatchField
(
    const calculatedProcessorFvPatchField<Type>& ptf
)
:
    coupledFvPatchField<Type>(ptf),
    procInterface_(ptf.procInterface_),
    sendBuf_(procInterface_.faceCells().size()),
    receiveBuf_(procInterface_.faceCells().size()),
    scalarSendBuf_(procInterface_.faceCells().size()),
    scalarReceiveBuf_(procInterface_.faceCells().size()),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)

{}


template<class Type>
Foam::calculatedProcessorFvPatchField<Type>::calculatedProcessorFvPatchField
(
    const calculatedProcessorFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    procInterface_(ptf.procInterface_),
    sendBuf_(procInterface_.faceCells().size()),
    receiveBuf_(procInterface_.faceCells().size()),
    scalarSendBuf_(procInterface_.faceCells().size()),
    scalarReceiveBuf_(procInterface_.faceCells().size()),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::calculatedProcessorFvPatchField<Type>::ready() const
{
    if
    (
        this->outstandingSendRequest_ >= 0
     && this->outstandingSendRequest_ < Pstream::nRequests()
    )
    {
        if (!UPstream::finishedRequest(this->outstandingSendRequest_))
        {
            return false;
        }
    }
    this->outstandingSendRequest_ = -1;

    if
    (
        this->outstandingRecvRequest_ >= 0
     && this->outstandingRecvRequest_ < Pstream::nRequests()
    )
    {
        if (!UPstream::finishedRequest(this->outstandingRecvRequest_))
        {
            return false;
        }
    }
    this->outstandingRecvRequest_ = -1;

    return true;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::calculatedProcessorFvPatchField<Type>::patchNeighbourField() const
{
    if (!this->ready())
    {
        FatalErrorInFunction
            << "On patch of size " << procInterface_.faceCells().size()
            << " between proc " << procInterface_.myProcNo()
            << " and " << procInterface_.neighbProcNo()
            << " outstanding request."
            << abort(FatalError);
    }
    return *this;
}


template<class Type>
void Foam::calculatedProcessorFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        if (!is_contiguous<Type>::value)
        {
            FatalErrorInFunction
                << "Invalid for non-contiguous data types"
                << abort(FatalError);
        }

        //this->patchInternalField(sendBuf_);
        // Bypass patchInternalField since uses fvPatch addressing
        {
            const Field<Type>& iF = this->internalField();
            const labelList& fc = procInterface_.faceCells();
            sendBuf_.setSize(fc.size());
            forAll(fc, i)
            {
                sendBuf_[i] = iF[fc[i]];
            }
        }

        // Receive straight into *this
        this->setSize(sendBuf_.size());
        outstandingRecvRequest_ = UPstream::nRequests();
        UIPstream::read
        (
            Pstream::commsTypes::nonBlocking,
            procInterface_.neighbProcNo(),
            this->data_bytes(),
            this->size_bytes(),
            procInterface_.tag(),
            procInterface_.comm()
        );

        outstandingSendRequest_ = UPstream::nRequests();
        UOPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            procInterface_.neighbProcNo(),
            sendBuf_.cdata_bytes(),
            sendBuf_.size_bytes(),
            procInterface_.tag(),
            procInterface_.comm()
        );
    }
}


template<class Type>
void Foam::calculatedProcessorFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        if
        (
            outstandingRecvRequest_ >= 0
         && outstandingRecvRequest_ < Pstream::nRequests()
        )
        {
            UPstream::waitRequest(outstandingRecvRequest_);
        }
        outstandingSendRequest_ = -1;
        outstandingRecvRequest_ = -1;
    }
}


template<class Type>
void Foam::calculatedProcessorFvPatchField<Type>::initInterfaceMatrixUpdate
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    // Bypass patchInternalField since uses fvPatch addressing

    const labelList& fc = lduAddr.patchAddr(patchId);

    scalarSendBuf_.setSize(fc.size());
    forAll(fc, i)
    {
        scalarSendBuf_[i] = psiInternal[fc[i]];
    }

    if (!this->ready())
    {
        FatalErrorInFunction
            << "On patch " //<< interface_.name()
            << " outstanding request."
            << abort(FatalError);
    }



    scalarReceiveBuf_.setSize(scalarSendBuf_.size());
    outstandingRecvRequest_ = UPstream::nRequests();

    UIPstream::read
    (
        Pstream::commsTypes::nonBlocking,
        procInterface_.neighbProcNo(),
        scalarReceiveBuf_.data_bytes(),
        scalarReceiveBuf_.size_bytes(),
        procInterface_.tag(),
        procInterface_.comm()
    );

    outstandingSendRequest_ = UPstream::nRequests();

    UOPstream::write
    (
        Pstream::commsTypes::nonBlocking,
        procInterface_.neighbProcNo(),
        scalarSendBuf_.cdata_bytes(),
        scalarSendBuf_.size_bytes(),
        procInterface_.tag(),
        procInterface_.comm()
    );

    const_cast<lduInterfaceField&>
    (
        static_cast<const lduInterfaceField&>(*this)
    ).updatedMatrix() = false;
}


template<class Type>
void Foam::calculatedProcessorFvPatchField<Type>::addToInternalField
(
    solveScalarField& result,
    const bool add,
    const scalarField& coeffs,
    const solveScalarField& vals
) const
{
    const labelUList& faceCells = this->procInterface_.faceCells();

    if (add)
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] += coeffs[elemI]*vals[elemI];
        }
    }
    else
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= coeffs[elemI]*vals[elemI];
        }
    }
}


template<class Type>
void Foam::calculatedProcessorFvPatchField<Type>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    if (this->updatedMatrix())
    {
        return;
    }


    if
    (
        outstandingRecvRequest_ >= 0
     && outstandingRecvRequest_ < Pstream::nRequests()
    )
    {
        UPstream::waitRequest(outstandingRecvRequest_);
    }
    // Recv finished so assume sending finished as well.
    outstandingSendRequest_ = -1;
    outstandingRecvRequest_ = -1;

    // Consume straight from scalarReceiveBuf_. Note use of our own
    // helper to avoid using fvPatch addressing
    addToInternalField(result, !add, coeffs, scalarReceiveBuf_);

    const_cast<lduInterfaceField&>
    (
        static_cast<const lduInterfaceField&>(*this)
    ).updatedMatrix() = true;
}


// ************************************************************************* //
