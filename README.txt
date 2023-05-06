ICSFoam is a library for Implicit Coupled Simulations in the finite-volume software OpenFOAM. A new density-based CFD solver is available for the simulation of high-speed flows.

The HiSA solver (J. Heyns, O. Oxtoby, 9th South African Conference on Computational and Applied Mechanics, SACAM 2014, https://gitlab.com/hisa/hisa) was uptaken as a starting point and was modified in order to make the code more compliant with the OF structure, to generalize the formulation to an arbitrary number of variables in the coupled system, and to make it suitable for state-of-art turbomachinery simulations. Support for fully-implicit Harmonic Balance Method simulations is implemented.

In order to compile and run the ICSFoam library and the applications, you'll need the v2112 ESI version of OpenFOAM. We cannot guarantee that the solver will compile and run on other versions. The files are distributed under the GNU GPL v3 license version or later.

There are a few OpenFOAM files that you need to replace to run the code in parallel. These are contained in orginalOFFiles/. In the repository main folder, after loading OF environment variables:

rm -r $FOAM_SRC/finiteVolume/fields/fvPatchFields/constraint/

cp -r originalOFFiles/constraintFvPatchFields/ $FOAM_SRC/finiteVolume/fields/fvPatchFields/constraint/

This will basically replace the constraintFvPatchFieds with the new ones (the only difference is just that a few private members have been made public). Then recompile the OF finiteVolume library with:

wmake libso /home/stefano/OpenFOAM/OpenFOAM-v2112/src/finiteVolume/

Finally, to compile the code, you need to run the Allwmake file in the repository main folder.

After that, it is strongly recommended to recompile Openfoam src and utilities to update dependecies (Allwmake inside the main OpenFoam directory).
For example, to decompose and reconstruct cases you need to recompile decomposePar and reconstructPar utilities.

If you use the code, please cite either of our papers: https://doi.org/10.1016/j.cpc.2023.108673, https://doi.org/10.3390/machines10040279

Many thanks to Dr. Nicola Casari and Dr. Mauro Carnevale for their help during the code implementation and test phases!

Stefano Oliani

