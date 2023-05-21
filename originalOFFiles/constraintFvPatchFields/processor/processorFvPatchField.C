/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "processorFvPatchField.H"
#include "processorFvPatch.H"
#include "demandDrivenData.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    procPatch_(refCast<const processorFvPatch>(p)),
    sendBuf_(0),
    receiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0)
{}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& f
)
:
    coupledFvPatchField<Type>(p, iF, f),
    procPatch_(refCast<const processorFvPatch>(p)),
    sendBuf_(0),
    receiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0)
{}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict, dict.found("value")),
    procPatch_(refCast<const processorFvPatch>(p, dict)),
    sendBuf_(0),
    receiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0)
{
    if (!isA<processorFvPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    // If the value is not supplied set to the internal field
    if (!dict.found("value"))
    {
        fvPatchField<Type>::operator=(this->patchInternalField());
    }
}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorFvPatch>(p)),
    sendBuf_(0),
    receiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0)
{
    if (!isA<processorFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }
    if (debug && !ptf.ready())
    {
        FatalErrorInFunction
            << "On patch " << procPatch_.name() << " outstanding request."
            << abort(FatalError);
    }
}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf
)
:
    processorLduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    procPatch_(refCast<const processorFvPatch>(ptf.patch())),
    sendBuf_(std::move(ptf.sendBuf_)),
    receiveBuf_(std::move(ptf.receiveBuf_)),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalarSendBuf_(std::move(ptf.scalarSendBuf_)),
    scalarReceiveBuf_(std::move(ptf.scalarReceiveBuf_))
{
    if (debug && !ptf.ready())
    {
        FatalErrorInFunction
            << "On patch " << procPatch_.name() << " outstanding request."
            << abort(FatalError);
    }
}


template<class Type>
Foam::processorFvPatchField<Type>::processorFvPatchField
(
    const processorFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorFvPatch>(ptf.patch())),
    sendBuf_(0),
    receiveBuf_(0),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1),
    scalarSendBuf_(0),
    scalarReceiveBuf_(0)
{
    if (debug && !ptf.ready())
    {
        FatalErrorInFunction
            << "On patch " << procPatch_.name() << " outstanding request."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::processorFvPatchField<Type>::patchNeighbourField() const
{
    if (debug && !this->ready())
    {
        FatalErrorInFunction
            << "On patch " << procPatch_.name()
            << " outstanding request."
            << abort(FatalError);
    }
    return *this;
}


template<class Type>
void Foam::processorFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        this->patchInternalField(sendBuf_);

        if
        (
            commsType == Pstream::commsTypes::nonBlocking
         && !Pstream::floatTransfer
        )
        {
            if (!is_contiguous<Type>::value)
            {
                FatalErrorInFunction
                    << "Invalid for non-contiguous data types"
                    << abort(FatalError);
            }

            // Receive straight into *this
            this->setSize(sendBuf_.size());
            outstandingRecvRequest_ = UPstream::nRequests();
            UIPstream::read
            (
                Pstream::commsTypes::nonBlocking,
                procPatch_.neighbProcNo(),
                this->data_bytes(),
                this->size_bytes(),
                procPatch_.tag(),
                procPatch_.comm()
            );

            outstandingSendRequest_ = UPstream::nRequests();
            UOPstream::write
            (
                Pstream::commsTypes::nonBlocking,
                procPatch_.neighbProcNo(),
                sendBuf_.cdata_bytes(),
                sendBuf_.size_bytes(),
                procPatch_.tag(),
                procPatch_.comm()
            );
        }
        else
        {
            procPatch_.compressedSend(commsType, sendBuf_);
        }
    }
}


template<class Type>
void Foam::processorFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        if
        (
            commsType == Pstream::commsTypes::nonBlocking
         && !Pstream::floatTransfer
        )
        {
            // Fast path. Received into *this

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
        else
        {
            procPatch_.compressedReceive<Type>(commsType, *this);
        }

        if (doTransform())
        {
            transform(*this, procPatch_.forwardT(), *this);
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::processorFvPatchField<Type>::snGrad
(
    const scalarField& deltaCoeffs
) const
{
    return deltaCoeffs*(*this - this->patchInternalField());
}


template<class Type>
void Foam::processorFvPatchField<Type>::initInterfaceMatrixUpdate
(
    solveScalarField&,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField&,
    const direction,
    const Pstream::commsTypes commsType
) const
{
    //this->patch().patchInternalField(psiInternal, scalarSendBuf_);

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    scalarSendBuf_.setSize(this->patch().size());
    forAll(scalarSendBuf_, facei)
    {
        scalarSendBuf_[facei] = psiInternal[faceCells[facei]];
    }

    if
    (
        commsType == Pstream::commsTypes::nonBlocking
     && !Pstream::floatTransfer
    )
    {
        // Fast path.
        if (debug && !this->ready())
        {
            FatalErrorInFunction
                << "On patch " << procPatch_.name()
                << " outstanding request."
                << abort(FatalError);
        }


        scalarReceiveBuf_.setSize(scalarSendBuf_.size());
        outstandingRecvRequest_ = UPstream::nRequests();
        UIPstream::read
        (
            Pstream::commsTypes::nonBlocking,
            procPatch_.neighbProcNo(),
            scalarReceiveBuf_.data_bytes(),
            scalarReceiveBuf_.size_bytes(),
            procPatch_.tag(),
            procPatch_.comm()
        );

        outstandingSendRequest_ = UPstream::nRequests();
        UOPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            procPatch_.neighbProcNo(),
            scalarSendBuf_.cdata_bytes(),
            scalarSendBuf_.size_bytes(),
            procPatch_.tag(),
            procPatch_.comm()
        );
    }
    else
    {
        procPatch_.compressedSend(commsType, scalarSendBuf_);
    }

    const_cast<processorFvPatchField<Type>&>(*this).updatedMatrix() = false;
}


template<class Type>
void Foam::processorFvPatchField<Type>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    if (this->updatedMatrix())
    {
        return;
    }

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    if
    (
        commsType == Pstream::commsTypes::nonBlocking
     && !Pstream::floatTransfer
    )
    {
        // Fast path.
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
        // Consume straight from scalarReceiveBuf_

        if (!std::is_arithmetic<Type>::value)
        {
            // Transform non-scalar data according to the transformation tensor
            transformCoupleField(scalarReceiveBuf_, cmpt);
        }

        // Multiply the field by coefficients and add into the result
        this->addToInternalField
        (
            result,
            !add,
            faceCells,
            coeffs,
            scalarReceiveBuf_
        );
    }
    else
    {
        solveScalarField pnf
        (
            procPatch_.compressedReceive<solveScalar>
            (
                commsType,
                this->size()
            )()
        );

        if (!std::is_arithmetic<Type>::value)
        {
            // Transform non-scalar data according to the transformation tensor
            transformCoupleField(pnf, cmpt);
        }

        // Multiply the field by coefficients and add into the result
        this->addToInternalField(result, !add, faceCells, coeffs, pnf);
    }

    const_cast<processorFvPatchField<Type>&>(*this).updatedMatrix() = true;
}


template<class Type>
void Foam::processorFvPatchField<Type>::initInterfaceMatrixUpdate
(
    Field<Type>&,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>& psiInternal,
    const scalarField&,
    const Pstream::commsTypes commsType
) const
{
    sendBuf_.setSize(this->patch().size());

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    forAll(sendBuf_, facei)
    {
        sendBuf_[facei] = psiInternal[faceCells[facei]];
    }

    if
    (
        commsType == Pstream::commsTypes::nonBlocking
     && !Pstream::floatTransfer
    )
    {
        // Fast path.
        if (debug && !this->ready())
        {
            FatalErrorInFunction
                << "On patch " << procPatch_.name()
                << " outstanding request."
                << abort(FatalError);
        }


        receiveBuf_.setSize(sendBuf_.size());
        outstandingRecvRequest_ = UPstream::nRequests();
        IPstream::read
        (
            Pstream::commsTypes::nonBlocking,
            procPatch_.neighbProcNo(),
            receiveBuf_.data_bytes(),
            receiveBuf_.size_bytes(),
            procPatch_.tag(),
            procPatch_.comm()
        );

        outstandingSendRequest_ = UPstream::nRequests();
        OPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            procPatch_.neighbProcNo(),
            sendBuf_.cdata_bytes(),
            sendBuf_.size_bytes(),
            procPatch_.tag(),
            procPatch_.comm()
        );
    }
    else
    {
        procPatch_.compressedSend(commsType, sendBuf_);
    }

    const_cast<processorFvPatchField<Type>&>(*this).updatedMatrix() = false;
}


template<class Type>
void Foam::processorFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>&,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{
    if (this->updatedMatrix())
    {
        return;
    }

    const labelUList& faceCells = lduAddr.patchAddr(patchId);

    if
    (
        commsType == Pstream::commsTypes::nonBlocking
     && !Pstream::floatTransfer
    )
    {
        // Fast path.
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

        // Consume straight from receiveBuf_

        // Transform according to the transformation tensor
        transformCoupleField(receiveBuf_);

        // Multiply the field by coefficients and add into the result
        this->addToInternalField(result, !add, faceCells, coeffs, receiveBuf_);
    }
    else
    {
        Field<Type> pnf
        (
            procPatch_.compressedReceive<Type>(commsType, this->size())()
        );

        // Transform according to the transformation tensor
        transformCoupleField(pnf);

        // Multiply the field by coefficients and add into the result
        this->addToInternalField(result, !add, faceCells, coeffs, pnf);
    }

    const_cast<processorFvPatchField<Type>&>(*this).updatedMatrix() = true;
}


template<class Type>
bool Foam::processorFvPatchField<Type>::ready() const
{
    if
    (
        outstandingSendRequest_ >= 0
     && outstandingSendRequest_ < Pstream::nRequests()
    )
    {
        bool finished = UPstream::finishedRequest(outstandingSendRequest_);
        if (!finished)
        {
            return false;
        }
    }
    outstandingSendRequest_ = -1;

    if
    (
        outstandingRecvRequest_ >= 0
     && outstandingRecvRequest_ < Pstream::nRequests()
    )
    {
        bool finished = UPstream::finishedRequest(outstandingRecvRequest_);
        if (!finished)
        {
            return false;
        }
    }
    outstandingRecvRequest_ = -1;

    return true;
}


// ************************************************************************* //
