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

Class
    Foam::overlapAMIFvsPatchField

Description
    Adapted from overlapGGI patch implemented in foam-extend by Hrvoje Jasak, 
    Fethi Tekin and Oliver Borm

SourceFiles
    overlapAMIFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef overlapAMIFvsPatchField_H
#define overlapAMIFvsPatchField_H

#include "coupledFvsPatchField.H"
#include "overlapAMIFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class overlapAMIFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class overlapAMIFvsPatchField
:
    public coupledFvsPatchField<Type>
{
    // Private data

        //- Local reference cast into the cyclic patch
        const overlapAMIFvPatch& overlapAMIPatch_;


public:

    //- Runtime type information
    TypeName(overlapAMIFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        overlapAMIFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        overlapAMIFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given overlapAMIFvsPatchField onto a new patch
        overlapAMIFvsPatchField
        (
            const overlapAMIFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        overlapAMIFvsPatchField
        (
            const overlapAMIFvsPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvsPatchField<Type>> clone() const
        {
            return tmp<fvsPatchField<Type>>
            (
                new overlapAMIFvsPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        overlapAMIFvsPatchField
        (
            const overlapAMIFvsPatchField<Type>&,
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
                new overlapAMIFvsPatchField<Type>(*this, iF)
            );
        }

    // Member functions

        // Access

            //- Return true if running parallel
            virtual bool coupled() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "overlapAMIFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
