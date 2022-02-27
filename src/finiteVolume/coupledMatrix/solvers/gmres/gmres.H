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

#ifndef gmres_H
#define gmres_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "residualsIO.H"
#include "fvCFD.H"
#include "coupledMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class gmresSolver Declaration
\*---------------------------------------------------------------------------*/

class gmres
:
    public coupledMatrix::solver
{
    // Number of Krylov-space vectors
    label nDirs_;

    protected:

        inline void givensRotation(const scalar& h, const scalar& beta, scalar& c, scalar& s) const;

    public:
        //- Runtime type information
        TypeName("GMRES");

        gmres
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

