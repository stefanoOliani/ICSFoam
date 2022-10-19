ICSFoam is a library for Implicit Coupled Simulations in the finite-volume software OpenFOAM. A new density-based CFD solver is available for the simulation of high-speed flows.

In order to compile and run the ICSFoam library and the applications, you'll need the v2006 version of OpenFOAM. We cannot guarantee that the solver will compile and run on other versions. The files are distributed under the GNU GPL v3 license version or later.

There are a few OpenFOAM files that you need to replace to run the code in parallel. These are contained in orginalOFFiles/. In the repository main folder, after loading OF environment variables:

rm -r $FOAM_SRC/finiteVolume/fields/fvPatchFields/constraint/

cp -r originalOFFiles/constraintFvPatchFields/ $FOAM_SRC/finiteVolume/fields/fvPatchFields/constraint/

This will basically replace the constraintFvPatchFieds with the new ones (the only difference is just that a few private members have been made public). Then recompile the OF finiteVolume library with:

wmake libso /home/stefano/OpenFOAM/OpenFOAM-v2006/src/finiteVolume/

Finally, to compile the code, you need to run the Allwmake file in the repository main folder.

After that, it is strongly recommended to recompile Openfoam src and utilities to update dependecies (Allwmake inside the main OpenFoam directory).
For example, to decompose and reconstruct cases you need to recompile decomposePar and reconstructPar utilities.


Many thanks to Dr. Nicola Casari and Dr. Mauro Carnevale for their help during the code implementation and test phases!

Stefano Oliani

