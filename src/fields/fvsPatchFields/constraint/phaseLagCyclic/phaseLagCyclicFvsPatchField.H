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


Class
    Foam::phaseLagCyclicFvsPatchField

Description
    Foam::phaseLagCyclicFvsPatchField

SourceFiles
    phaseLagCyclicFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef phaseLagCyclicFvsPatchField_H
#define phaseLagCyclicFvsPatchField_H

#include "cyclicFvsPatchField.H"
#include "phaseLagCyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseLagCyclicFvsPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class phaseLagCyclicFvsPatchField
:
    public cyclicFvsPatchField<Type>
{
    // Private data

        //- Local reference cast into the phaseLagCyclic patch
        const phaseLagCyclicFvPatch& phaseLagCyclicPatch_;


public:

    //- Runtime type information
    TypeName("phaseLagCyclic");


    // Constructors

        //- Construct from patch and internal field
        phaseLagCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        phaseLagCyclicFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given phaseLagCyclicFvsPatchField onto a new patch
        phaseLagCyclicFvsPatchField
        (
            const phaseLagCyclicFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        phaseLagCyclicFvsPatchField
        (
            const phaseLagCyclicFvsPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvsPatchField<Type>> clone() const
        {
            return tmp<fvsPatchField<Type>>
            (
                new phaseLagCyclicFvsPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        phaseLagCyclicFvsPatchField
        (
            const phaseLagCyclicFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvsPatchField<Type>> clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return tmp<fvsPatchField<Type>>
            (
                new phaseLagCyclicFvsPatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "phaseLagCyclicFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
