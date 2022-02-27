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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::pseudotimeControl::nCorrOuter() const
{
    return nCorrOuter_;
}


inline Foam::residualsIO& Foam::pseudotimeControl::residuals()
{
    return residuals_;
}


inline bool Foam::pseudotimeControl::storeInitialResiduals() const
{
    // Start from second outer iteration
    return (corr_ == 2);
}


inline bool Foam::pseudotimeControl::firstIter() const
{
    return corr_ == 1;
}


inline bool Foam::pseudotimeControl::finalIter() const
{
    return converged_ || (corr_ == nCorrOuter_);
}


inline bool Foam::pseudotimeControl::turbCorr() const
{
    return !turbOnFinalIterOnly_ || finalIter();
}


inline void Foam::pseudotimeControl::setResidual(const residualsIO& res)
{
    residuals_ = res;
}


// ************************************************************************* //
