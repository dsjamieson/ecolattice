
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2020			
	 *
	 *	 main function of the simulation program. instantiates an
	 *	 Ecolattice object and calls 'run.' 
	 *
	 ***********************************************************/

#include "ecolattice.hpp"

int main(int argc, char* argv[]) {
	if (argc != 2) { 
		fprintf(stderr, "Usage: ecolattice input_filename\n");
		exit(0);
	}
	Ecolattice sim(argv[1], 0);
	sim.run();
	return(0);
}

