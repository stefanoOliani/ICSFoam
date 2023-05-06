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
    Foam::overlapAMIGAMGInterfaceField

Description
    Adapted from overlapGGI patch implemented in foam-extend by Hrvoje Jasak, 
    Fethi Tekin and Oliver Borm

SourceFiles
    overlapAMIGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef overlapAMIGAMGInterfaceField_H
#define overlapAMIGAMGInterfaceField_H

#include "GAMGInterfaceField.H"
#include "overlapAMIGAMGInterface.H"
#include "overlapAMILduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class overlapAMIGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class overlapAMIGAMGInterfaceField
:
    public GAMGInterfaceField,
    virtual public overlapAMILduInterfaceField
{
    // Private data

        //- Local reference cast into the cyclic interface
        const overlapAMIGAMGInterface& overlapAMIInterface_;

        //- Is the transform required
        bool doTransform_;

        //- Rank of component for transformation
        int rank_;


    // Private Member Functions

        //- No copy construct
        overlapAMIGAMGInterfaceField
        (
            const overlapAMIGAMGInterfaceField&
        ) = delete;

        //- No copy assignment
        void operator=(const overlapAMIGAMGInterfaceField&) = delete;


public:

    //- Runtime type information
    TypeName("overlapAMI");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        overlapAMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );

        //- Construct from GAMG interface and fine level interface field
        overlapAMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const bool doTransform,
            const int rank
        );


    //- Destructor
    virtual ~overlapAMIGAMGInterfaceField();


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return overlapAMIInterface_.size();
            }


        // Interface matrix update

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                solveScalarField& result,
                const bool add,
				const lduAddressing& lduAddr,
				const label patchId,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;


        //- Cyclic interface functions

            //- Does the interface field perform the transformation
            virtual bool doTransform() const
            {
                return doTransform_;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return overlapAMIInterface_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return overlapAMIInterface_.reverseT();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return rank_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
