
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
#include "site_stepper.hpp"

void Simulation::run(void) {
	int start_time;
	int persistence;
	std::chrono::steady_clock::time_point t_start, t_end;
	std::chrono::duration<double, std::milli> duration;
	
	start_time = getContinueTime() + 1;  // start time can be 1 or any other number (less than duration) for restarts (repeats)
	// save initial state at time 0, including parameters, competition matrices, and initial positions of individuals
	if (start_time == 1) {
		saveCompetition();
		saveLattice(0);
	}

	#pragma omp parallel
	{
		#pragma omp single
		{
			t_start = std::chrono::steady_clock::now();
			fprintf(stdout, "Setting up thread RNGs\n");
		}
		// each thread gets a copy of the same RNG
		// random count initialized to zero
		//int my_id = omp_get_thread_num();
		std::mt19937 this_random_generator;
		seedGenerator(this_random_generator);
		unsigned long long this_random_count =  getRandomCount();
		discardRandom(getMaxRandomCount(), this_random_count, this_random_generator);
		this_random_count = getRandomCount();
		SiteStepper stepper(*this, this_random_generator, this_random_count);
		int i, j;
		#pragma omp single
		{
			t_end = std::chrono::steady_clock::now();
			duration = (t_end - t_start) / 1000.;
			fprintf(stdout, "Done setting up RNGs in %.4f seconds\n", duration.count());
			fprintf(stdout, "Starting at time step %d\n", start_time);
		}
		// run simulation
		for (int time_step = start_time; time_step < max_time_step + 1; time_step++) {
			#pragma omp single
			t_start = std::chrono::steady_clock::now();

			// distribute sites in lattice among threads
			#pragma omp for schedule(monotonic:static)
			for (int index = 0; index < lattice_size * lattice_size; index++) {
		
				i = index / lattice_size;
				j = index - i * lattice_size;
				// RNG discards values so that the simulations are repeatable
				discardRandom(4 * lattice_size * lattice_size * ((unsigned long long) time_step - 1) + 4 * (j + lattice_size * i) - this_random_count, 
		  						  this_random_count, this_random_generator);
				// update single site in 'next_lattice' and 'next_dispersal_lattice'
				//sim.updateSingleSite(i, j, this_random_count, this_random_generator);
				stepper.updateSingleSite(i, j);
			}

			#pragma omp single
			persistence = getPersistence();  // get the number of species currently coexisting in the lattice
			// if persistence is below the minimum allowed persistence, new random seeds are drawn and the simulation is reinitialized
			// the new initial state of the simulation is saved
			if (persistence < min_persistence) {
				#pragma omp single
				{
					fprintf(stdout, "\n\nToo many extinctions, reinitializing\n\n");
					setRandomSeeds();
					reinitializeSimulation(time_step);
					saveCompetition();
					saveLattice(0);
					persistence = 0;
					time_step = 0;
				}
				// pass new seeds to each thread
				seedGenerator(this_random_generator);
				this_random_count = 0;
				discardRandom(getMaxRandomCount(), this_random_count, this_random_generator);
				this_random_count = getRandomCount();
			}
			else {
				// species and seed locations in 'next_lattice' and 'next_dispersal_lattice' are converted to locations in 'lattice' and 'dispersal_lattice'
				// t + 1 becomes t for the next time step
				nextToThis();
				// write species locations to file
				#pragma omp single nowait
				saveLattice(time_step);
				saveDispersal(time_step);
				#pragma omp single
				{
					t_end = std::chrono::steady_clock::now();
					duration = (t_end - t_start) / 1000.;
					fprintf(stdout, "Done step %d of %d in %.4f seconds\n", time_step, max_time_step, duration.count());
				}
			}
		}
	}
	return;
}


