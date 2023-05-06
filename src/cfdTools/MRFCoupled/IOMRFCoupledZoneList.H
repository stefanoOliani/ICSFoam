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
    Foam::IOMRFCoupledZoneList

Description
    List of MRF zones with IO functionality.  MRF zones are specified by a list
    of dictionary entries, e.g.

    \verbatim
    zone1
    {
        cellZone    rotor1;
        active      yes;
        ...
    }

    zone2
    {
        cellZone    rotor2;
        active      yes;
        ...
    }
    \endverbatim

SourceFiles
    IOMRFCoupledZoneList.C

\*---------------------------------------------------------------------------*/

#ifndef IOMRFCoupledZoneList_H
#define IOMRFCoupledZoneList_H

#include "IOdictionary.H"
#include "MRFCoupledZoneList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class IOMRFCoupledZoneList Declaration
\*---------------------------------------------------------------------------*/

class IOMRFCoupledZoneList
:
    public IOdictionary,
    public MRFCoupledZoneList
{
private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- No copy construct
        IOMRFCoupledZoneList(const IOMRFCoupledZoneList&) = delete;

        //- No copy assignment
        void operator=(const IOMRFCoupledZoneList&) = delete;


public:

    // Constructors

        //- Construct from mesh
        IOMRFCoupledZoneList(const fvMesh& mesh);


        //- Destructor
        virtual ~IOMRFCoupledZoneList() = default;


    // Member Functions

        //- Read dictionary
        virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
