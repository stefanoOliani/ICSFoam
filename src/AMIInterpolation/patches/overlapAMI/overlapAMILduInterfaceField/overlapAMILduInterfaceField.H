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
    Foam::overlapAMILduInterfaceField

Description
    Adapted from overlapGGI patch implemented in foam-extend by Hrvoje Jasak, 
    Fethi Tekin and Oliver Borm

SourceFiles
    overlapAMILduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef overlapAMILduInterfaceField_H
#define overlapAMILduInterfaceField_H

#include "primitiveFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class overlapAMILduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class overlapAMILduInterfaceField
{

public:

    //- Runtime type information
    TypeName("overlapAMILduInterfaceField");


    // Constructors

        //- Construct null
        overlapAMILduInterfaceField()
        {}


    //- Destructor
    virtual ~overlapAMILduInterfaceField();


    // Member Functions

        // Access

            //- Is the transform required
            virtual bool doTransform() const = 0;

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const = 0;

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const = 0;

            //- Return rank of component for transform
            virtual int rank() const = 0;


        //- Transform given patch field
        template<class Type>
        void transformCoupleField(Field<Type>& f) const;

        //- Transform given patch internal field
        void transformCoupleField
        (
            solveScalarField& psiInternal,
            const direction cmpt
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tensorField.H"

template<class Type>
void Foam::overlapAMILduInterfaceField::transformCoupleField
(
    Field<Type>& f
) const
{
    if (doTransform())
    {
        if (forwardT().size() == 1)
        {
            transform(f, forwardT()[0], f);
        }
        else
        {
            transform(f, forwardT(), f);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
