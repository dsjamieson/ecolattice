# Options are MPI, OMP, (both parallel) or SERIES (not parallel)
TYPE=SERIES
# Installation directory
INSTALLDIR=bin
 
ifeq ($(TYPE), MPI)
	# MPI c++ compiler
	CXX=mpic++ 
	CXXFLAGS=-std=c++11
	SOURCE=$(wildcard src/mpi/*.cpp)
else ifeq ($(TYPE), OMP)
	# c++ compiler
	CXX=g++ 
	CXXFLAGS=-O3 -msse3 -fpermissive -fopenmp -std=c++11
	SOURCE=$(wildcard src/omp/*.cpp)
else ifeq ($(TYPE), SERIES)
	# c++ compiler
	CXX=g++ 
	CXXFLAGS=-std=c++11
	SOURCE=$(wildcard src/series/*.cpp)
endif

EXEC=ecolattice

all: ecolattice

ecolattice: $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(SOURCE)

install: $(EXEC)
	mkdir -p $(INSTALLDIR)
	cp $(EXEC) $(INSTALLDIR)

uninstall: 
	rm -f $(INSTALLDIR)/$(EXEC)

clean:
	rm -f $(EXEC)


