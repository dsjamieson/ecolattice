# Options are PARALLEL (with MPI) or SERIES (without MPI)
TYPE = PARALLEL
# Installation directory
INSTALLDIR = bin 
 
ifeq ($(TYPE), PARALLEL)
	# MPI c++ compiler
	CXX = mpic++ 
	SOURCE = $(wildcard src/parallel/*.cpp)
else ifeq ($(TYPE), SERIES)
	# c++ compiler
	CXX = g++ 
	SOURCE = $(wildcard src/series/*.cpp)
endif

CXXFLAGS = -std=c++11
EXEC = ecolattice

all: ecolat

ecolat: $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(SOURCE)

install: $(EXEC)
	mkdir -p $(INSTALLDIR)
	cp $(EXEC) $(INSTALLDIR)

uninstall: 
	rm -f $(INSTALLDIR)/$(EXEC)

clean:
	rm -f $(EXEC)

