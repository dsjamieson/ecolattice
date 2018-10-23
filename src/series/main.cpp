
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2018			
	 *
	 *	 this is the main body of the simulation program if
	 *	 you do not want to run in parallel and do not have
	 *	 MPI installed.
	 *	 
	 *	
	 *
	 ***********************************************************/

#include "simulation.h"

int main(int argc, char* argv[]) {

	if (argc != 2) { 
		fprintf(stderr, "Usage: ecolattice input_filename\n");
		exit(0);
	}

	int i, j, box_size, time_step, max_time_step;
	
	// instantiating class Simulation with object sim
	// requires file name with list of parameters
	Simulation sim(argv[1], 0);
	
	// save initial state at time 0, including parameters, competition matrices, and initial positions of individuals
	sim.saveCompetition();
	sim.saveProperties();
	sim.saveBox(0);

	box_size = sim.getBoxSize();  // number of cells in one dimension of the box
	max_time_step = sim.getMaxTimeStep();  // duration of the simulation

	// run simulation
	for (time_step = 1; time_step < max_time_step + 1; time_step++) {
		for (i = 0; i < box_size; i++) {
			for (j = 0;j < box_size; j++) {
				// RNG discards values so that the simulations are repeatable
				sim.discardRandom(((int64_t) (4 * box_size * box_size * (time_step - 1) + 4 * (j + box_size * i))) - sim.getRandomCount());
				// update single site in 'next_box' and 'next_dispersal_box'
				sim.updateSingleSite(i, j);
			}
		}
		
		// species and seed locations in 'next_box' and 'next_dispersal_box' are converted to locations in 'box' and 'dispersal_box'
		// t + 1 becomes t for the next time step
		sim.nextToThis();
		// write species locations to file
		sim.saveBox(time_step);
	}
	return(0);
}


