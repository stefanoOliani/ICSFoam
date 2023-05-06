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
    Foam::phaseLagAMIGAMGInterfaceField

Description
    GAMG agglomerated cyclic interface field.

SourceFiles
    phaseLagAMIGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef phaseLagAMIGAMGInterfaceField_H
#define phaseLagAMIGAMGInterfaceField_H

#include "GAMGInterfaceField.H"
#include "phaseLagAMIGAMGInterface.H"
#include "phaseLagAMILduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class phaseLagAMIGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class phaseLagAMIGAMGInterfaceField
:
    public GAMGInterfaceField,
    virtual public phaseLagAMILduInterfaceField
{
    // Private data

        //- Local reference cast into the cyclic interface
        const phaseLagAMIGAMGInterface& phaseLagAMIInterface_;

        //- Is the transform required
        bool doTransform_;

        //- Rank of component for transformation
        int rank_;


    // Private Member Functions

        //- No copy construct
        phaseLagAMIGAMGInterfaceField
        (
            const phaseLagAMIGAMGInterfaceField&
        ) = delete;

        //- No copy assignment
        void operator=(const phaseLagAMIGAMGInterfaceField&) = delete;


public:

    //- Runtime type information
    TypeName("phaseLagAMI");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        phaseLagAMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );

        //- Construct from GAMG interface and fine level interface field
        phaseLagAMIGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const bool doTransform,
            const int rank
        );


    //- Destructor
    virtual ~phaseLagAMIGAMGInterfaceField();


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return phaseLagAMIInterface_.size();
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
                return phaseLagAMIInterface_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return phaseLagAMIInterface_.reverseT();
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
