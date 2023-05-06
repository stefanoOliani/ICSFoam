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
    Foam::IOHBZoneList

Description
    List of HB zones with IO functionality.  HB zones are specified by a list
    of dictionary entries, e.g.

    \verbatim
    zone1
    {
        cellZone    rotor1;
        active      yes;
        frequencies 2(3.14 6.28)
        ...
    }

    zone2
    {
        cellZone    rotor2;
        active      yes;
        frequencies 2(1 3.3)
        ...
    }
    \endverbatim

SourceFiles
    IOHBZoneList.C

\*---------------------------------------------------------------------------*/

#ifndef IOHBZoneList_H
#define IOHBZoneList_H

#include "HBZoneList.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class IOHBZoneList Declaration
\*---------------------------------------------------------------------------*/

class IOHBZoneList
:
    public IOdictionary,
    public HBZoneList
{
private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- No copy construct
        IOHBZoneList(const IOHBZoneList&) = delete;

        //- No copy assignment
        void operator=(const IOHBZoneList&) = delete;


public:

    // Constructors

        //- Construct from mesh
        IOHBZoneList(const fvMesh& mesh);


        //- Destructor
        virtual ~IOHBZoneList() = default;


    // Member Functions

        //- Read dictionary
        virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
