/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa

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
    lusgs

Description
    LU-SGS preconditioning as per Luo, Baum, Lohner, JCP 146 664-690 (1998).

Authors


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
