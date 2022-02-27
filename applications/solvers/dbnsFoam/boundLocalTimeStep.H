   if (localTimesteppingBounding)
    {
        volScalarField rhoMin =
            localTimesteppingLowerBound*scalarVarsPrevIter[0];

        volScalarField eMin =
            localTimesteppingLowerBound*
            (scalarVarsPrevIter[1]/scalarVarsPrevIter[0]
            - 0.5*magSqr(vectorVarsPrevIter[0]/scalarVarsPrevIter[0]));

        volScalarField eTemp(rhoE/rho - 0.5*magSqr(rhoU/rho));

        volScalarField factor
        (
            IOobject
            (
                "factor",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("factor", dimless, 1.0),
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(mesh.owner(), facei)
        {
            label own = mesh.owner()[facei];
            label nei = mesh.neighbour()[facei];

            if
            (
                (rho[own] < rhoMin[own]) 
             || (eTemp[own] < eMin[own]) 
             || (eTemp[own] < SMALL)
            )
            {
                factor[own] = min(0.5, factor[own]);
                factor[nei] = min(0.75, factor[nei]);
            }
            if
            (
                (rho[nei] < rhoMin[nei]) 
             || (eTemp[nei] < eMin[nei]) 
             || (eTemp[nei] < SMALL)
            )
            {
                factor[nei] = min(0.5, factor[nei]);
                factor[own] = min(0.75, factor[own]);
            }
        }

        forAll(mesh.boundary(), patchi)
        {
            if (mesh.boundary()[patchi].coupled())
            {
                scalarField rhoNei =
                    rho.boundaryField()[patchi].patchNeighbourField();

                scalarField rhoMinNei =
                    rhoMin.boundaryField()[patchi].patchNeighbourField();

                scalarField eTempNei =
                    eTemp.boundaryField()[patchi].patchNeighbourField();

                scalarField eMinNei =
                    eMin.boundaryField()[patchi].patchNeighbourField();

                const labelUList& fc = mesh.boundary()[patchi].faceCells();

                forAll(fc, bfacei)
                {
                    if
                    (
                        (rho[fc[bfacei]] < rhoMin[fc[bfacei]])
                     || (eTemp[fc[bfacei]] < eMin[fc[bfacei]])
                     || (eTemp[fc[bfacei]] < SMALL)
                    )
                    {
                        factor[fc[bfacei]] = min(0.5, factor[fc[bfacei]]);
                    }
                    if
                    (
                        (rhoNei[bfacei] < rhoMinNei[bfacei])
                     || (eTempNei[bfacei] < eMinNei[bfacei])
                     || (eTempNei[bfacei] < SMALL)
                    )
                    {
                        factor[fc[bfacei]] = min(0.75, factor[fc[bfacei]]);
                    }
                }
            }
        }

        pseudoCoField *= factor;
    }

