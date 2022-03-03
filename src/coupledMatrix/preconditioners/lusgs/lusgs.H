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

#ifndef lusgs_H
#define lusgs_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#include "coupledMatrix.H"
#include "fvCFD.H"
#include "blockFvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class lusgs Declaration
\*---------------------------------------------------------------------------*/

class lusgs
:
    public coupledMatrix::preconditioner
{

    protected:

        //- Protected data
            autoPtr<scalarField> rDiagCoeff_;

        //- Private member functions

            template <class Type1, class Type2>
            void forwardSweep
            (
                const label celli,
                Field<Type1>& result,
                const blockFvMatrix<Type1, typename outerProduct<Type1, Type2>::type>& matrix,
                const Type2& delta
            ) const;

            template <class Type1, class Type2>
            void reverseSweep
            (
                const label celli,
                Field<Type1>& result,
                const blockFvMatrix<Type1, typename outerProduct<Type1, Type2>::type>& matrix,
                const Type2& delta
            ) const;

            void divideByDiagonal
            (
                const label& celli,
                scalarList& dS,
                List<vector>& dV,
                const PtrList<scalarField>& sVec,
                const PtrList<vectorField>& vVec
            ) const;

    public:
        //- Runtime type information
        TypeName("LUSGS");

        lusgs
        (
            const coupledMatrix::solver& sol,
            const dictionary& solverControlsUnused
        );

        //- Precondition a vector
        virtual void precondition
        (
            PtrList<scalarField>& sVec,
            PtrList<vectorField>& vVec
        ) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
