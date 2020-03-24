
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2018			
	 *
	 *	 this is the main body of the simulation program to run
	 *	 with OMP. ten threads are initialized, and each 
	 *	 thread works through the lattice in parallel. 
	 *	 may only be used with one node.
	 *
	 ***********************************************************/

#include "simulation.hpp"

int main(int argc, char* argv[]) {
	if (argc != 2) { 
		fprintf(stderr, "Usage: ecolattice input_filename\n");
		exit(0);
	}
	Simulation sim(argv[1], 0);
	sim.run();
	return(0);
}

