
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

	int i, j, lattice_size, time_step, start_time, max_time_step;
	int persistence, min_persistence;
	
	// instantiating class Simulation with object sim
	// requires file name with list of parameters
	Simulation sim(argv[1], 0);

	lattice_size = sim.getLatticeSize();  // number of cells in one dimension of the lattice
	start_time = sim.getRestartTime() + 1;
	max_time_step = sim.getMaxTimeStep();  // duration of the simulation
	min_persistence = sim.getMinPersistence();

	// save initial state at time 0, including parameters, competition matrices, and initial positions of individuals
	if (start_time == 1) {
		sim.saveCompetition();
		sim.saveLattice(0);
	}

	// run simulation
	for (time_step = start_time; time_step < max_time_step + 1; time_step++) {

		for (i = 0; i < lattice_size; i++) {
			for (j = 0;j < lattice_size; j++) {
				// RNG discards values so that the simulations are repeatable
				sim.discardRandom(((unsigned long long) (4 * lattice_size * lattice_size * (time_step - 1) + 4 * (j + lattice_size * i))) - sim.getRandomCount());
				// update single site in 'next_lattice' and 'next_dispersal_lattice'
				sim.updateSingleSite(i, j);
			}
		}
		
		persistence = sim.getPersistence();  // get the number of species currently coexisting in the lattice
		
		// if persistence is below the minimum allowed persistence, new random seeds are drawn and the simulation is reinitialized. the new initial state of the simulation is save
		if ( persistence < min_persistence ) {
			fprintf(stdout, "Too many extinctions, reinitializing\n");
			sim.setRandomSeeds();
			sim.reinitializeSimulation(time_step);
			sim.saveCompetition();
			sim.saveLattice(0);
			persistence = 0;
			time_step = 0;
			continue;
		}
		
		// species and seed locations in 'next_lattice' and 'next_dispersal_lattice' are converted to locations in 'lattice' and 'dispersal_lattice'
		// t + 1 becomes t for the next time step
		sim.nextToThis();
		// write species locations to file
		sim.saveLattice(time_step);
		sim.saveDispersal(time_step);

		fprintf(stdout, "done step %d\n", time_step);
	}
	return(0);
}


