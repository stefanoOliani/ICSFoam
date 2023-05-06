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
    smoothSolverCoupled

Description


Authors

\*---------------------------------------------------------------------------*/

#ifndef smoothSolverCoupled_H
#define smoothSolverCoupled_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "residualsIO.H"
#include "fvCFD.H"
#include "coupledMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class smoothSolverCoupledSolver Declaration
\*---------------------------------------------------------------------------*/

class smoothSolverCoupled
:
    public coupledMatrix::solver
{
protected:

    // Protected data

        //- Number of sweeps before the evaluation of residual
        label nSweeps_;

        //- Read the control parameters from the controlDict_
        virtual void readControls();

    public:
        //- Runtime type information
        TypeName("smoothSolverCoupled");

        smoothSolverCoupled
        (
            const dictionary& dict,
            const coupledMatrix& matrix
        );

        // Main solve functions:

        //- Returns solver performance
        virtual residualsIO solve
        (
            PtrList<volScalarField>& sW, PtrList<volVectorField>& vW,
            const PtrList<scalarField>& sSource, const PtrList<vectorField>& vSource
        ) const;

        //- Returns solver performance
        virtual residualsIO solveDelta
        (
            PtrList<volScalarField>& sW, PtrList<volVectorField>& vW,
            const PtrList<scalarField>& sSource, const PtrList<vectorField>& vSource
        ) const;

        //- Returns solver performance
        virtual residualsIO solveDelta
        (
            PtrList<volScalarField>& sW, PtrList<volVectorField>& vW,
            const PtrList<scalarField>& sSource, const PtrList<vectorField>& vSource,
			PtrList<volScalarField>& dsW, PtrList<volVectorField>& dvW
        ) const;


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

