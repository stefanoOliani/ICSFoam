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
    Foam::overlapAMILduInterface

Description
    Adapted from overlapGGI patch implemented in foam-extend by Hrvoje Jasak, 
    Fethi Tekin and Oliver Borm

SourceFiles
    overlapAMILduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef overlapAMILduInterface_H
#define overlapAMILduInterface_H

#include "primitiveFieldsFwd.H"
#include "AMIPatchToPatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class overlapAMILduInterface Declaration
\*---------------------------------------------------------------------------*/

class overlapAMILduInterface
{

public:

    //- Runtime type information
    TypeName("overlapAMILduInterface");


    // Constructors

        //- Construct null
        overlapAMILduInterface() = default;


    //- Destructor
    virtual ~overlapAMILduInterface() = default;


    // Member Functions

        // Access

            //- Return neighbour
            virtual label neighbPatchID() const = 0;

            virtual bool owner() const = 0;

            //- Return processor number
            virtual const overlapAMILduInterface& neighbPatch() const = 0;

            virtual const AMIPatchToPatchInterpolation& AMI() const = 0;

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const = 0;

            //- Return face reverse transformation tensor
            virtual const tensorField& reverseT() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
