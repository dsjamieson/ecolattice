# ecolattice

a C++ program to conduct stochastic spatial lattice simulations of plant communities. ecolattice can run in series on a single processor, or can run in parallel using OpenMP or MPI.

## authors

Drew S. Jamieson (drew.s.jamieson@gmail.com) and Nicole L. Kinlock (nlkinlock@gmail.com)

## system requirements

C++11 compliant compiler (GCC is used in `make`)

**for MPI version of ecolattice** MPICH 3.1.4

**for OMP version of ecolattice** OpenMP 4.5

## installation

install from github

```sh
git clone https://github.com/dsjamieson/ecolattice
```

edit `Makefile` as needed (`$TYPE`) for single-processor version (`SERIES`) or parallel versions (`OMP` or `MPI`). then, to compile the program, install the program in desired directory (may be changed in `Makefile`, `$INSTALLDIR`)

```sh
make install clean
```

## usage

input files must follow the structure of the sample input file, which can be found in `ecolattice/example/parameters.dat`. create directory where output will be written (OutfileDir) before running ecolattice.

**running in series or with OMP**

```sh
ecolattice/bin/ecolattice inputfile
```

**running with MPI**

```sh
mpirun -np <number processors> ecolattice/bin/ecolattice inputfile
```
