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

#include "blockFvMatrix.H"
#include "LduInterfaceField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class sourceType, class blockType>
template<class psiType>
void Foam::blockFvMatrix<sourceType, blockType>::initProcessorInterfaces
(
	const Field<psiType>& psiif,
	Field<psiType>& result,
	const LduInterfaceFieldPtrsList<psiType>& interfaces
) const
{
	scalarField interfaceDummyCoeffs(1);

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(interfaces, interfacei)
        {
            if (interfaces.set(interfacei))
            {
                interfaces[interfacei].initInterfaceMatrixUpdate
                (
                    result,
                    true,
                    psiif,
					interfaceDummyCoeffs,
                    Pstream::defaultCommsType
                );
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfacei=patchSchedule.size()/2;
            interfacei<interfaces.size();
            interfacei++
        )
        {
            if (interfaces.set(interfacei))
            {
                interfaces[interfacei].initInterfaceMatrixUpdate
                (
                    result,
                    true,
                    psiif,
					interfaceDummyCoeffs,
                    Pstream::commsTypes::blocking
                );
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}


template<class sourceType, class blockType>
template<class psiType>
void Foam::blockFvMatrix<sourceType, blockType>::updateProcessorInterfaces
(
	Field<sourceType>& result,
	const Pstream::commsTypes commsType,
	const processorFvPatchField<psiType>& procField,
	const fvMesh& mesh,
	const label patchi
) const
{
    if (procField.updatedMatrix())
    {
        return;
    }

    const processorFvPatch& procPatch = refCast<const processorFvPatch>(procField.patch());

//    if(commsType == Pstream::commsTypes::scheduled)
//    {
//    	Info<<"scheduled "<<endl;
//    }
//
//    if(commsType == Pstream::commsTypes::blocking)
//    {
//    	Info<<"blocking "<<endl;
//    }
//
//    if(commsType == Pstream::commsTypes::nonBlocking)
//    {
//    	Info<<"non blocking "<<endl;
//    }

    if
    (
        commsType == Pstream::commsTypes::nonBlocking
     && !Pstream::floatTransfer
    )
    {
        // Fast path.
        if
        (
            procField.outstandingRecvRequest() >= 0
         && procField.outstandingRecvRequest() < Pstream::nRequests()
        )
        {
            UPstream::waitRequest(procField.outstandingRecvRequest());
        }
        // Recv finished so assume sending finished as well.
        procField.outstandingSendRequest() = -1;
        procField.outstandingRecvRequest() = -1;

        // Consume straight from receiveBuf_

        // Transform according to the transformation tensor
        procField.transformCoupleField(procField.receiveBuf());

        const labelUList& faceCells = mesh.boundary()[patchi].faceCells();

        // Diagonal matrix does not have interfaces defined
        if (this->interfacesUpper().size())
        {
    	    const Field<blockType>& pdSByS = this->interfacesUpper()[patchi];

    	    forAll(procField.receiveBuf(), facei)
    	    {
    		    result[faceCells[facei]] += dot(pdSByS[facei],procField.receiveBuf()[facei]);
    	    }
        }
    }
    else
    {
        Field<psiType> pnf
        (
            procPatch.compressedReceive<psiType>(commsType, procField.size())()
        );

        // Transform according to the transformation tensor
        procField.transformCoupleField(pnf);

        const labelUList& faceCells = mesh.boundary()[patchi].faceCells();

        // Diagonal matrix does not have interfaces defined
        if (this->interfacesUpper().size())
        {
    	    const Field<blockType>& pdSByS = this->interfacesUpper()[patchi];

    	    forAll(pnf, facei)
    	    {
    		    result[faceCells[facei]] += dot(pdSByS[facei],pnf[facei]);
    	    }
        }
    }

    const_cast<processorFvPatchField<psiType>&>(procField).updatedMatrix() = true;
}
// ************************************************************************* //
