# Installation directory
INSTALLDIR=bin
 
# c++ compiler
CXX=g++
OPT += -DOMP
CXXFLAGS= -Wall -O2 -msse3 -fpermissive -fopenmp -std=c++17 $(OPT)
SOURCE=$(wildcard src/*.cpp)

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


