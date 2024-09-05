To calculate absorption, circular dichroism (CD), emission and Circularly Polarized Luminescence (CPL) spectra using the exciton model, just follow the steps below:

Prepare the supramolecular geometry file in the xyz format. A template (positions.dat) is provided.
Prepare a file with the xyz components of electric transition dipole moments. A template (electric-dipoles.dat) is provided.
Prepare a file with the xyz components of magnetic transition dipole moments. A template (magnetic-dipoles.dat) is provided.

Compile absorption.f90 if you are interested in computing absorption and CD, or emission.f90 if you are interested in computing emission and CPL, possibly using the Intel Fortran compiler (MKL routines are required):

ifort absorption.f90 -o absorption.e -mkl

ifort emission.f90 -o emission.e -mkl


Run the executable(s):

./absorption.e

./emission.e

The programs will ask for additional inputs, which can be inserted from terminal.

The outputs are two data files which can be plotted to show the spectra. For absorption.f90 the output files are absorption.dat (contains the UV-Vis absorption spectrum) and CD.dat (contains the UV-Vis CD spectrum); for emission.f90 the output files are emission.dat (contains the emission spectrum) and CPL.dat (contains CPL spectrum)