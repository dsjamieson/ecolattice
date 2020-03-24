
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

#include "simulation.h"

int main(int argc, char* argv[]) {

	if (argc != 2) { 
		fprintf(stderr, "Usage: ecolattice input_filename\n");
		exit(0);
	}

	// determine number of processors and set thread IDs
	int nproc = 1;
	nproc = 1;
	omp_set_num_threads(nproc);
	thread_local int my_id;
	#pragma omp parallel
	{
		my_id = omp_get_thread_num();
	}
	fprintf(stdout, "Ecolattice running in parallel (OMP) with %d threads\n", nproc);

	int lattice_size, time_step, start_time, max_time_step;
	int persistence, min_persistence;
	
	// instantiating class Simulation with object sim
	// requires file name with list of parameters
	Simulation sim(argv[1], 0);
	// each thread gets a copy of the same RNG
	// random count initialized to zero
	thread_local std::mt19937 this_random_generator;
	thread_local unsigned long long this_random_count = 0;
	#pragma omp parallel 
	{
		sim.seedGenerator(this_random_generator);
		sim.discardRandom(sim.getMaxRandomCount(), this_random_count, this_random_generator);
		this_random_count = sim.getRandomCount();
	}

	lattice_size = sim.getLatticeSize();  // number of cells in one dimension of the lattice
	start_time = sim.getRestartTime() + 1;  // start time can be 1 or any other number (less than duration) for restarts (repeats)
	max_time_step = sim.getMaxTimeStep();  // duration of the simulation
	min_persistence = sim.getMinPersistence();  // minimum species richness before reinitializing simulation

	// save initial state at time 0, including parameters, competition matrices, and initial positions of individuals
	if (start_time == 1) {
		sim.saveCompetition();
		sim.saveLattice(0);
	}

	fprintf(stdout, "Starting at time step %d\n", start_time);

	// run simulation
	for (time_step = start_time; time_step < max_time_step + 1; time_step++) {

		std::chrono::steady_clock::time_point t_start = std::chrono::steady_clock::now();

		// distribute sites in lattice among threads
		#pragma omp parallel for schedule(monotonic:static) 
		for (int index = 0; index < lattice_size * lattice_size; index++) {
		
			int i = index / lattice_size;
			int j = index - i * lattice_size;

			// RNG discards values so that the simulations are repeatable
			sim.discardRandom(4 * lattice_size * lattice_size * ((unsigned long long) time_step - 1) + 4 * (j + lattice_size * i) - this_random_count, 
								this_random_count, this_random_generator);

			// update single site in 'next_lattice' and 'next_dispersal_lattice'
			sim.updateSingleSite(i, j, this_random_count, this_random_generator);
		}

		persistence = sim.getPersistence();  // get the number of species currently coexisting in the lattice
		
		// if persistence is below the minimum allowed persistence, new random seeds are drawn and the simulation is reinitialized
		// the new initial state of the simulation is saved
		if (persistence < min_persistence) {
			fprintf(stdout, "\n\nToo many extinctions, reinitializing\n\n");
			sim.setRandomSeeds();
			sim.reinitializeSimulation(time_step);
			sim.saveCompetition();
			sim.saveLattice(0);
			persistence = 0;
			time_step = 0;
			// pass new seeds to each thread
			#pragma omp parallel 
			{
				sim.seedGenerator(this_random_generator);
				this_random_count = 0;
				sim.discardRandom(sim.getMaxRandomCount(), this_random_count, this_random_generator);
				this_random_count = sim.getRandomCount();
			}
			continue;
		}
		
		// species and seed locations in 'next_lattice' and 'next_dispersal_lattice' are converted to locations in 'lattice' and 'dispersal_lattice'
		// t + 1 becomes t for the next time step
		sim.nextToThis();
		// write species locations to file
		sim.saveLattice(time_step);
		sim.saveDispersal(time_step);

		std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
		std::chrono::duration<double, std::milli> duration = (t_end - t_start) / 1000.;
		fprintf(stdout, "Done step %d of %d in %.4f seconds\n", time_step, max_time_step, duration);

	}
	return(0);
}


