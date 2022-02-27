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

#ifndef characteristicBase_H
#define characteristicBase_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class characteristicBase Declaration
\*---------------------------------------------------------------------------*/

class characteristicBase
{
protected:
    // Protected data

        //- Name of velocity field, default = "U"
        word UName_;

        //- Name of velocity field, default = "p"
        word pName_;

        //- Name of temperature field, default = "T"
        word TName_;

        //- Name of mass flux field, default = "phi"
        word phiName_;

        //- Name of density field, default = "rho"
        word rhoName_;

        //- External (specified) velocity
        vector URef_;

        //- External (specified) pressure
        scalar pRef_;

        //- External (specified) temperature
        scalar TRef_;


public:

    //- Runtime type information
    TypeName("characteristicBase");


    // Constructors

        //- Construct from patch and internal field
        characteristicBase
        (
            const fvPatch&
        );

        //- Construct from patch, internal field and dictionary
        characteristicBase
        (
            const fvPatch&,
            const dictionary&
        );

        //- Construct by mapping given characteristicBase
        //  onto a new patch
        characteristicBase
        (
            const characteristicBase&,
            const fvPatch&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        characteristicBase
        (
            const characteristicBase&
        );

    //- Destructor
    virtual ~characteristicBase() 
    {};


    // Member functions

        // Access

            //- Return the external velocity
            const vector& URef() const
            {
                return URef_;
            }

            //- Return reference to the external velocity to allow adjustment
            vector& URef()
            {
                return URef_;
            }

            //- Return the external pressure
            scalar pRef() const
            {
                return pRef_;
            }

            //- Return reference to the external pressure to allow adjustment
            scalar& pRef()
            {
                return pRef_;
            }

            //- Return the temperature at infinity
            scalar TRef() const
            {
                return TRef_;
            }

            //- Return reference to the external temperature
            //  to allow adjustment
            scalar& TRef()
            {
                return TRef_;
            }


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
