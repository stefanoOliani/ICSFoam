/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa

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


#include "blockFvMatrix.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class sourceType, class blockType>
template<class Type2>
void Foam::blockFvMatrix<sourceType,blockType>::addToInternalField
(
    const labelUList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    if (addr.size() != pf.size())
    {
        FatalErrorInFunction
            << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    forAll(addr, facei)
    {
        intf[addr[facei]] += pf[facei];
    }
}

template<class sourceType, class blockType>
template<class Type2>
void Foam::blockFvMatrix<sourceType,blockType>::addToInternalField
(
    const labelUList& addr,
    const tmp<Field<Type2>>& tpf,
    Field<Type2>& intf
) const
{
    addToInternalField(addr, tpf(), intf);
    tpf.clear();
}


template<class sourceType, class blockType>
template<class Type2>
void Foam::blockFvMatrix<sourceType,blockType>::subtractFromInternalField
(
    const labelUList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    if (addr.size() != pf.size())
    {
        FatalErrorInFunction
            << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    forAll(addr, facei)
    {
        intf[addr[facei]] -= pf[facei];
    }
}


template<class sourceType, class blockType>
template<class Type2>
void Foam::blockFvMatrix<sourceType,blockType>::subtractFromInternalField
(
    const labelUList& addr,
    const tmp<Field<Type2>>& tpf,
    Field<Type2>& intf
) const
{
    subtractFromInternalField(addr, tpf(), intf);
    tpf.clear();
}


template<class sourceType, class blockType>
void Foam::blockFvMatrix<sourceType,blockType>::addBoundaryDiag
(
	fvMatrix<sourceType> matrix,
    scalarField& diag,
    const direction solveCmpt
) const
{
    forAll(matrix.internalCoeffs(), patchi)
    {
        addToInternalField
        (
            this->lduAddr().patchAddr(patchi),
            matrix.internalCoeffs()[patchi].component(solveCmpt),
            diag
        );
    }
}


template<class sourceType, class blockType>
void Foam::blockFvMatrix<sourceType,blockType>::addBoundarySource
(
	fvMatrix<sourceType> matrix,
    Field<sourceType>& source,
    const bool couples
) const
{
    forAll(matrix.psi().boundaryField(), patchi)
    {
        const fvPatchField<sourceType>& ptf = matrix.psi().boundaryField()[patchi];
        const Field<sourceType>& pbc = matrix.boundaryCoeffs()[patchi];

        if (!ptf.coupled())
        {
            addToInternalField(this->lduAddr().patchAddr(patchi), pbc, source);
        }
        else if (couples)
        {
            const tmp<Field<sourceType>> tpnf = ptf.patchNeighbourField();
            const Field<sourceType>& pnf = tpnf();

            const labelUList& addr = this->lduAddr().patchAddr(patchi);

            forAll(addr, facei)
            {
                source[addr[facei]] += cmptMultiply(pbc[facei], pnf[facei]);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class sourceType, class blockType>
Foam::blockFvMatrix<sourceType,blockType>::blockFvMatrix
(
    const lduMesh& mesh
)
:
    refCount(),
    LduMatrix<sourceType,blockType,blockType>(mesh)
{}

template<class sourceType, class blockType>
Foam::blockFvMatrix<sourceType,blockType>::blockFvMatrix(const blockFvMatrix<sourceType,blockType>& bfm)
:
    refCount(),
    LduMatrix<sourceType,blockType,blockType>(bfm)
{
    // Perform assignment of interfacesUpper and interfacesLower, which doesn't happen in LduMatrix
    this->interfacesUpper().resize(bfm.interfacesUpper().size());
    this->interfacesLower().resize(bfm.interfacesUpper().size());

    forAll(this->interfacesUpper(), i)
    {
        this->interfacesUpper().set(i, new Field<blockType>(bfm.interfacesUpper()[i]));
    }
    forAll(this->interfacesLower(), i)
    {
        this->interfacesLower().set(i, new Field<blockType>(bfm.interfacesLower()[i]));
    }
}

template<class sourceType, class blockType>
Foam::tmp<Foam::blockFvMatrix<sourceType,blockType>> Foam::blockFvMatrix<sourceType,blockType>::clone() const
{
    return tmp<blockFvMatrix<sourceType,blockType>>
    (
        new blockFvMatrix<sourceType,blockType>(*this)
    );
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class sourceType, class blockType>
Foam::blockFvMatrix<sourceType,blockType>::~blockFvMatrix()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class sourceType, class blockType>
void Foam::blockFvMatrix<sourceType,blockType>::insertBlock
(
	const GeometricField<blockType, fvsPatchField, surfaceMesh>& leftField,
	const GeometricField<blockType, fvsPatchField, surfaceMesh>& rightField
)
{
	const fvMesh& mesh = leftField.mesh();

    const labelUList& owner = mesh.owner();

    blockFvMatrix<sourceType, blockType> mx(mesh);

    Field<blockType>& upp = mx.upper();
    Field<blockType>& low = mx.lower();
    Field<blockType>& diag = mx.diag();

    forAll(owner, facei)
    {
       	scalar magSfi = mesh.magSf()[facei];

	    upp[facei] = 0.5 * magSfi * rightField[facei];
	    low[facei] = -0.5 * magSfi * leftField[facei];
    }

    mx.negSumDiag();

    mx.interfacesUpper().resize(mesh.boundary().size());
    mx.interfacesLower().resize(mesh.boundary().size());

    forAll(mesh.boundary(), patchi)
    {
       	const scalarField& magSfp = mesh.magSf().boundaryField()[patchi];

    	const Field<blockType>& leftFp = leftField.boundaryField()[patchi];
    	const Field<blockType>& rightFp = rightField.boundaryField()[patchi];

	    mx.interfacesUpper().set(patchi, 0.5*magSfp * rightFp);
	    mx.interfacesLower().set(patchi, -0.5*magSfp * leftFp);

	   // Don't include physical boundaries because those are dealt with separately
	   if (mesh.boundary()[patchi].coupled())
	   {
		   const labelUList& bOwner = mesh.boundary()[patchi].faceCells();

           Field<blockType>& low = mx.interfacesLower()[patchi];

		   forAll(bOwner, facei)
		   {
			   const label own = bOwner[facei];

			   diag[own] -= low[facei];
		   }
	   }
    }

    operator+=(mx);
}


template<class sourceType, class blockType>
void Foam::blockFvMatrix<sourceType,blockType>::insertDissipationBlock
(
	const GeometricField<blockType, fvsPatchField, surfaceMesh>& dissField
)
{
	const fvMesh& mesh = dissField.mesh();

    const labelUList& owner = mesh.owner();

    blockFvMatrix<sourceType, blockType> mx(mesh);

    Field<blockType>& upp = mx.upper();
    Field<blockType>& low = mx.lower();
    Field<blockType>& diag = mx.diag();

    forAll(owner, facei)
    {
       	scalar magSfi = mesh.magSf()[facei];

	    upp[facei] = 0.5 * magSfi * dissField[facei];
	    low[facei] = 0.5 * magSfi * dissField[facei];
    }

    mx.negSumDiag();

    mx.interfacesUpper().resize(mesh.boundary().size());
    mx.interfacesLower().resize(mesh.boundary().size());

    forAll(mesh.boundary(), patchi)
    {
       	const scalarField& magSfp = mesh.magSf().boundaryField()[patchi];

    	const Field<blockType>& dissFp = dissField.boundaryField()[patchi];

	    mx.interfacesUpper().set(patchi, 0.5*magSfp * dissFp);
	    mx.interfacesLower().set(patchi, 0.5*magSfp * dissFp);

	    //Include physical boundaries or not?
//	   if (mesh.boundary()[patchi].coupled())
//	   {
		   const labelUList& bOwner = mesh.boundary()[patchi].faceCells();

		   Field<blockType>& low = mx.interfacesLower()[patchi];

		   forAll(bOwner, facei)
		   {
			   const label own = bOwner[facei];

			   diag[own] -= low[facei];
		   }
//	   }
    }

    operator-=(mx);
}


template<class sourceType, class blockType>
template<class psiType>
void Foam::blockFvMatrix<sourceType,blockType>::Amul
(
    Field<sourceType>& Apsi,
	const GeometricField<psiType, fvPatchField, volMesh>& psi,
	const fvMesh& mesh
) const
{
	scalarField interfaceDummyCoeffs(1);

    sourceType* __restrict__ ApsiPtr = Apsi.begin();

    const psiType* const __restrict__ psiPtr = psi.primitiveField().begin();

    const blockType* const __restrict__ diagPtr = this->diag().begin();

    const label* const __restrict__ uPtr = this->lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = this->lduAddr().lowerAddr().begin();

    LduInterfaceFieldPtrsList<psiType> interfaces =
                psi.boundaryField().interfaces();

    Field<psiType> dummyResult(1);

    const label startRequest = Pstream::nRequests();

    // Initialise the update of interfaced interfaces
    initProcessorInterfaces
    (
		psi.primitiveField(),
        dummyResult,
		interfaces
    );

    const label nCells = this->diag().size();
    for (label cell=0; cell<nCells; cell++)
    {
        ApsiPtr[cell] = dot(diagPtr[cell], psiPtr[cell]);
    }

    if (this->hasLower() || this->hasUpper())
    {
        const blockType* const __restrict__ upperPtr = this->upper().begin();
        const blockType* const __restrict__ lowerPtr = this->lower().begin();

        const label nFaces = this->upper().size();
        for (label face=0; face<nFaces; face++)
        {
            ApsiPtr[uPtr[face]] += dot(lowerPtr[face], psiPtr[lPtr[face]]);
            ApsiPtr[lPtr[face]] += dot(upperPtr[face], psiPtr[uPtr[face]]);
        }
    }

    if
	(
		Pstream::defaultCommsType == Pstream::commsTypes::blocking
	)
	{
		forAll(mesh.boundary(), patchi)
		{
		   if (isA<processorFvPatch>(mesh.boundary()[patchi]))
		   {
			   const processorFvPatchField<psiType>& procField =
					   (refCast<const processorFvPatchField<psiType>>(psi.boundaryField()[patchi]));

			   updateProcessorInterfaces
			   (
					Apsi,
					Pstream::defaultCommsType,
					procField,
					mesh,
					patchi
			   );
		   }
		   else if (mesh.boundary()[patchi].coupled())
		   {
			   const labelUList& faceCells = mesh.boundary()[patchi].faceCells();

			   // Diagonal matrix does not have interfaces defined
			   if (this->interfacesUpper().size())
			   {
				   Field<psiType> pdS = psi.boundaryField()[patchi].patchNeighbourField();

				   const Field<blockType>& pdSByS = this->interfacesUpper()[patchi];

				   forAll(pdS, facei)
				   {
					   Apsi[faceCells[facei]] += dot(pdSByS[facei],pdS[facei]);
				   }
			   }
		   }
		}
	}
    else if
	(
		Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
	)
	{
    	// Try and consume interfaces as they become available
		bool allUpdated = false;

        // Block for everything
		if (Pstream::parRun())
		{
			if (allUpdated)
			{
				// All received. Just remove all outstanding requests
				UPstream::resetRequests(startRequest);
			}
			else
			{
				// Block for all requests and remove storage
				UPstream::waitRequests(startRequest);
			}
		}

		// Consume
		forAll(interfaces, interfacei)
		{
			if
			(
				interfaces.set(interfacei)
			&& !interfaces[interfacei].updatedMatrix()
			)
			{
			   if (isA<processorFvPatch>(mesh.boundary()[interfacei]))
			   {
				   const processorFvPatchField<psiType>& procField =
						   (refCast<const processorFvPatchField<psiType>>(psi.boundaryField()[interfacei]));

				   updateProcessorInterfaces
				   (
						Apsi,
						Pstream::defaultCommsType,
						procField,
						mesh,
						interfacei
				   );
			   }
			   else if (mesh.boundary()[interfacei].coupled())
			   {
				   const labelUList& faceCells = mesh.boundary()[interfacei].faceCells();

				   // Diagonal matrix does not have interfaces defined
				   if (this->interfacesUpper().size())
				   {
					   Field<psiType> pdS = psi.boundaryField()[interfacei].patchNeighbourField();

					   const Field<blockType>& pdSByS = this->interfacesUpper()[interfacei];

					   forAll(pdS, facei)
					   {
						   Apsi[faceCells[facei]] += dot(pdSByS[facei],pdS[facei]);
					   }
				   }
			   }
		   }
		}
	}
	else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
	{
		const lduSchedule& patchSchedule = this->patchSchedule();

		// Loop over all the "normal" interfaces relating to standard patches
		forAll(patchSchedule, i)
		{
			label interfacei = patchSchedule[i].patch;

			   if (mesh.boundary()[interfacei].coupled())
			   {
				   if (patchSchedule[i].init)
				   {
					   interfaces[interfacei].initInterfaceMatrixUpdate
					   (
						   dummyResult,
						   true,
						   this->lduAddr(),
						   interfacei,
						   psi.primitiveField(),
						   interfaceDummyCoeffs,
						   Pstream::defaultCommsType
					   );
				   }
				   else
				   {
					   if(isA<processorFvPatch>(mesh.boundary()[interfacei]))
					   {
						   const processorFvPatchField<psiType>& procField =
								   (refCast<const processorFvPatchField<psiType>>(psi.boundaryField()[interfacei]));

						   updateProcessorInterfaces
						   (
								Apsi,
								Pstream::defaultCommsType,
								procField,
								mesh,
								interfacei
						   );
					   }
					   else
					   {
						   const labelUList& faceCells = mesh.boundary()[interfacei].faceCells();

							   // Diagonal matrix does not have interfaces defined
							   if (this->interfacesUpper().size())
							   {
								   Field<psiType> pdS = psi.boundaryField()[interfacei].patchNeighbourField();

								   const Field<blockType>& pdSByS = this->interfacesUpper()[interfacei];

								   forAll(pdS, facei)
								   {
									   Apsi[faceCells[facei]] += dot(pdSByS[facei],pdS[facei]);
								   }
							   }
					   }
				   }
			   }
		}

		// Loop over the "global" patches are on the list of interfaces but
		// beyond the end of the schedule which only handles "normal" patches
		for
		(
			label interfacei=patchSchedule.size()/2;
			interfacei<interfaces.size();
			interfacei++
		)
		{

			if(isA<processorFvPatch>(mesh.boundary()[interfacei]))
			{
			   const processorFvPatchField<psiType>& procField =
					   (refCast<const processorFvPatchField<psiType>>(psi.boundaryField()[interfacei]));

			   updateProcessorInterfaces
			   (
					Apsi,
					Pstream::defaultCommsType,
					procField,
					mesh,
					interfacei
			   );
			}
			else if (mesh.boundary()[interfacei].coupled())
			{
			   const labelUList& faceCells = mesh.boundary()[interfacei].faceCells();

			   // Diagonal matrix does not have interfaces defined
			   if (this->interfacesUpper().size())
			   {
					Field<psiType> pdS = psi.boundaryField()[interfacei].patchNeighbourField();

					const Field<blockType>& pdSByS = this->interfacesUpper()[interfacei];

					forAll(pdS, facei)
					{
					   Apsi[faceCells[facei]] += dot(pdSByS[facei],pdS[facei]);
					}
				}
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
void Foam::blockFvMatrix<sourceType,blockType>::AmulNoDiag
(
    Field<sourceType>& Apsi,
	const GeometricField<psiType, fvPatchField, volMesh>& psi,
	const fvMesh& mesh
) const
{
	Apsi = Zero;
	scalarField interfaceDummyCoeffs(1);

    sourceType* __restrict__ ApsiPtr = Apsi.begin();

    const psiType* const __restrict__ psiPtr = psi.primitiveField().begin();

    const label* const __restrict__ uPtr = this->lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = this->lduAddr().lowerAddr().begin();

    LduInterfaceFieldPtrsList<psiType> interfaces =
                psi.boundaryField().interfaces();

    Field<psiType> dummyResult(1);

    const label startRequest = Pstream::nRequests();

    // Initialise the update of interfaced interfaces
    initProcessorInterfaces
    (
		psi.primitiveField(),
        dummyResult,
		interfaces
    );

    if (this->hasLower() || this->hasUpper())
    {
        const blockType* const __restrict__ upperPtr = this->upper().begin();
        const blockType* const __restrict__ lowerPtr = this->lower().begin();

        const label nFaces = this->upper().size();

        for (label face=0; face<nFaces; face++)
        {
            ApsiPtr[uPtr[face]] += dot(lowerPtr[face], psiPtr[lPtr[face]]);
            ApsiPtr[lPtr[face]] += dot(upperPtr[face], psiPtr[uPtr[face]]);
        }
    }

    if
	(
		Pstream::defaultCommsType == Pstream::commsTypes::blocking
	)
	{
		forAll(mesh.boundary(), patchi)
		{
		   if (isA<processorFvPatch>(mesh.boundary()[patchi]))
		   {
			   const processorFvPatchField<psiType>& procField =
					   (refCast<const processorFvPatchField<psiType>>(psi.boundaryField()[patchi]));

			   updateProcessorInterfaces
			   (
					Apsi,
					Pstream::defaultCommsType,
					procField,
					mesh,
					patchi
			   );
		   }
		   else if (mesh.boundary()[patchi].coupled())
		   {
			   const labelUList& faceCells = mesh.boundary()[patchi].faceCells();

			   // Diagonal matrix does not have interfaces defined
			   if (this->interfacesUpper().size())
			   {
				   Field<psiType> pdS = psi.boundaryField()[patchi].patchNeighbourField();

				   const Field<blockType>& pdSByS = this->interfacesUpper()[patchi];

				   forAll(pdS, facei)
				   {
					   Apsi[faceCells[facei]] += dot(pdSByS[facei],pdS[facei]);
				   }
			   }
		   }
		}
	}
    else if
	(
		Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
	)
	{
    	// Try and consume interfaces as they become available
		bool allUpdated = false;

        // Block for everything
		if (Pstream::parRun())
		{
			if (allUpdated)
			{
				// All received. Just remove all outstanding requests
				UPstream::resetRequests(startRequest);
			}
			else
			{
				// Block for all requests and remove storage
				UPstream::waitRequests(startRequest);
			}
		}

		// Consume
		forAll(interfaces, interfacei)
		{
			if
			(
				interfaces.set(interfacei)
			&& !interfaces[interfacei].updatedMatrix()
			)
			{
			   if (isA<processorFvPatch>(mesh.boundary()[interfacei]))
			   {
				   const processorFvPatchField<psiType>& procField =
						   (refCast<const processorFvPatchField<psiType>>(psi.boundaryField()[interfacei]));

				   updateProcessorInterfaces
				   (
						Apsi,
						Pstream::defaultCommsType,
						procField,
						mesh,
						interfacei
				   );
			   }
			   else if (mesh.boundary()[interfacei].coupled())
			   {
				   const labelUList& faceCells = mesh.boundary()[interfacei].faceCells();

				   // Diagonal matrix does not have interfaces defined
				   if (this->interfacesUpper().size())
				   {
					   Field<psiType> pdS = psi.boundaryField()[interfacei].patchNeighbourField();

					   const Field<blockType>& pdSByS = this->interfacesUpper()[interfacei];

					   forAll(pdS, facei)
					   {
						   Apsi[faceCells[facei]] += dot(pdSByS[facei],pdS[facei]);
					   }
				   }
			   }
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
void Foam::blockFvMatrix<sourceType,blockType>::insertEquation
(
    fvMatrix<sourceType>& matrix
)
{}


template<class sourceType, class blockType>
void Foam::blockFvMatrix<sourceType, blockType>::operator=(const blockFvMatrix<sourceType, blockType>& blfm)
{
    LduMatrix<sourceType,blockType,blockType>::operator=(blfm);

    // Perform assignment of interfacesUpper and interfacesLower, which doesn't happen in LduMatrix
    this->interfacesUpper().resize(blfm.interfacesUpper().size());
    this->interfacesLower().resize(blfm.interfacesUpper().size());

    forAll(this->interfacesUpper(), i)
    {
        this->interfacesUpper().set(i, new Field<blockType>(blfm.interfacesUpper()[i]));
    }
    forAll(this->interfacesLower(), i)
    {
        this->interfacesLower().set(i, new Field<blockType>(blfm.interfacesLower()[i]));
    }
}

template<class sourceType, class blockType>
void Foam::blockFvMatrix<sourceType, blockType>::operator=(const tmp<blockFvMatrix<sourceType, blockType> >& tblfm)
{
    operator=(tblfm());
    tblfm.clear();
}

template<class sourceType, class blockType>
Foam::blockFvMatrix<sourceType, blockType>& Foam::blockFvMatrix<sourceType, blockType>::operator+=
(
	const blockFvMatrix<sourceType,
	blockType>& B
)
{
    // If B has extra interfacesUpper or interfacesLower not present here,
    // we have to expand to accommodate them before calling LduMatrix::operator+=
    const label nBupper = B.interfacesUpper().size();
    const label nupper = this->interfacesUpper().size();
    if (nBupper > nupper)
    {
        this->interfacesUpper().resize(nBupper);
        for(label i = nupper; i < nBupper; i++)
        {
            this->interfacesUpper().set(i,new Field<blockType>(B.interfacesUpper()[i].size(), pTraits<blockType>::zero));
        }
    }
    const label nBlower = B.interfacesLower().size();
    const label nlower = this->interfacesLower().size();
    if (nBlower > nlower)
    {
        this->interfacesLower().resize(nBlower);
        for(label i = nlower; i < nBlower; i++)
        {
            this->interfacesLower().set(i, new Field<blockType>(B.interfacesLower()[i].size(), pTraits<blockType>::zero));
        }
    }

    LduMatrix<sourceType, blockType, blockType>::operator+=(B);

    return *this;
}

template<class sourceType, class blockType>
Foam::blockFvMatrix<sourceType, blockType>& Foam::blockFvMatrix<sourceType, blockType>::operator+=
(
	const tmp<blockFvMatrix<sourceType, blockType> >& tB
)
{
    operator+=(tB());
    tB.clear();
    return *this;
}

template<class sourceType, class blockType>
Foam::blockFvMatrix<sourceType, blockType>& Foam::blockFvMatrix<sourceType, blockType>::operator-=
(
	const blockFvMatrix<sourceType, blockType>& B
)
{
    // If B has extra interfacesUpper or interfacesLower not present here,
    // we have to expand to accommodate them before calling LduMatrix::operator+=
    const label nBupper = B.interfacesUpper().size();
    const label nupper = this->interfacesUpper().size();
    if (nBupper > nupper)
    {
        this->interfacesUpper().resize(nBupper);
        for(label i = nupper; i < nBupper; i++)
        {
            this->interfacesUpper().set(i, new Field<blockType>(B.interfacesUpper()[i].size(), pTraits<blockType>::zero));
        }
    }
    const label nBlower = B.interfacesLower().size();
    const label nlower = this->interfacesLower().size();
    if (nBlower > nlower)
    {
        this->interfacesLower().resize(nBlower);
        for(label i = nlower; i < nBlower; i++)
        {
            this->interfacesLower().set(i, new Field<blockType>(B.interfacesLower()[i].size(), pTraits<blockType>::zero));
        }
    }

    LduMatrix<sourceType, blockType, blockType>::operator-=(B);
    return *this;
}

template<class sourceType, class blockType>
Foam::blockFvMatrix<sourceType, blockType>& Foam::blockFvMatrix<sourceType, blockType>::operator-=
(
	const tmp<blockFvMatrix<sourceType,blockType> >& tB
)
{
    operator-=(tB());
    tB.clear();
    return *this;
}


template<class sourceType, class blockType>
Foam::tmp<Foam::blockFvMatrix<sourceType, blockType> > Foam::operator+
(
    const blockFvMatrix<sourceType, blockType>& A,
    const blockFvMatrix<sourceType, blockType>& B
)
{
    tmp<blockFvMatrix<sourceType, blockType> > tC(new blockFvMatrix<sourceType, blockType>(A));
    tC() += B;
    return tC;
}


template<class sourceType, class blockType>
Foam::tmp<Foam::blockFvMatrix<sourceType, blockType> > Foam::operator+
(
    const tmp<blockFvMatrix<sourceType, blockType> >& tA,
    const blockFvMatrix<sourceType, blockType>& B
)
{
    tmp<blockFvMatrix<sourceType, blockType> > tC(tA.ptr());
    tC() += B;
    return tC;
}


template<class sourceType, class blockType>
Foam::tmp<Foam::blockFvMatrix<sourceType, blockType> > Foam::operator+
(
    const blockFvMatrix<sourceType, blockType>& A,
    const tmp<blockFvMatrix<sourceType, blockType> >& tB
)
{
    tmp<blockFvMatrix<sourceType, blockType> > tC(tB.ptr());
    tC() += A;
    return tC;
}


template<class sourceType, class blockType>
Foam::tmp<Foam::blockFvMatrix<sourceType, blockType> > Foam::operator+
(
    const tmp<blockFvMatrix<sourceType, blockType> >& tA,
    const tmp<blockFvMatrix<sourceType, blockType> >& tB
)
{
    tmp<blockFvMatrix<sourceType, blockType> > tC(tA.ptr());
    tC.ref() += tB();
    tB.clear();
    return tC;
}


template<class sourceType, class blockType>
Foam::tmp<Foam::blockFvMatrix<sourceType, blockType> > Foam::operator-
(
    const blockFvMatrix<sourceType, blockType>& A,
    const blockFvMatrix<sourceType, blockType>& B
)
{
    tmp<blockFvMatrix<sourceType, blockType> > tC(new blockFvMatrix<sourceType, blockType>(A));
    tC() -= B;
    return tC;
}


template<class sourceType, class blockType>
Foam::tmp<Foam::blockFvMatrix<sourceType, blockType> > Foam::operator-
(
    const tmp<blockFvMatrix<sourceType, blockType> >& tA,
    const blockFvMatrix<sourceType, blockType>& B
)
{
    tmp<blockFvMatrix<sourceType, blockType> > tC(tA.ptr());
    tC() -= B;
    return tC;
}


template<class sourceType, class blockType>
Foam::tmp<Foam::blockFvMatrix<sourceType, blockType> > Foam::operator-
(
    const blockFvMatrix<sourceType, blockType>& A,
    const tmp<blockFvMatrix<sourceType, blockType> >& tB
)
{
    tmp<blockFvMatrix<sourceType, blockType> > tC(tB.ptr());
    tC() -= A;
    return tC;
}


template<class sourceType, class blockType>
Foam::tmp<Foam::blockFvMatrix<sourceType, blockType> > Foam::operator-
(
    const tmp<blockFvMatrix<sourceType, blockType> >& tA,
    const tmp<blockFvMatrix<sourceType, blockType> >& tB
)
{
    tmp<blockFvMatrix<sourceType, blockType> > tC(tA.ptr());
    tC() -= tB();
    tB.clear();
    return tC;
}


//- Multiplication by a uniform value (a primitive)

template<class Type2>
Foam::tmp<Foam::blockFvMatrix<Foam::vector, Type2> > Foam::operator*(const tmp<blockFvMatrix<scalar, scalar> >& tA, const Type2& B)
{
    tmp<blockFvMatrix<vector, Type2> > tC(new blockFvMatrix<vector, Type2>(tA->mesh()));

    if (tA->hasDiag())
    {
        tC->diag() = tA->diag()*B;
    }

    if (tA->hasUpper())
    {
        tC->upper() = tA->upper()*B;
    }

    if (tA->hasLower())
    {
        tC->lower() = tA->lower()*B;
    }

    tC->interfacesUpper() = tA->interfacesUpper()*B;
    tC->interfacesLower() = tA->interfacesLower()*B;

    tA.clear();
    return tC;
}

// ************************************************************************* //
