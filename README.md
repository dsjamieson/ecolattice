# ecolattice

a C++ program to conduct stochastic spatial lattice simulations of plant communities. ecolattice can run in series on a single processor, or can run in parallel using OpenMP.

## authors

Drew S. Jamieson (drew.s.jamieson@gmail.com) and Nicole L. Kinlock (nlkinlock@gmail.com)

## system requirements

C++11 compliant compiler (GCC is used in `make`)

**for OMP version of ecolattice** OpenMP 4.5

## installation

install from github

```sh
git clone https://github.com/dsjamieson/ecolattice
```

edit `Makefile` as needed by including or commenting out the `-DOMP` flag for parallel version or single-processor version, respectively. then, to compile the program, install in desired directory (may be changed in `Makefile`, `INSTALLDIR`)

```sh
make install clean
```

## usage

input files must follow the structure of the sample input file, which can be found in `ecolattice/example/parameters.dat`. create directory where output will be written (OutfileDir) before running ecolattice.

```sh
ecolattice/bin/ecolattice inputfile
```


