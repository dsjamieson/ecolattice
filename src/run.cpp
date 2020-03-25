
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2020			
	 *
	 *	driver of the Ecolattice object.  
	 *
	 ***********************************************************/

#include "ecolattice.hpp"
#include "site_stepper.hpp"

void Ecolattice::run(void) {
	int start_time;
	int persistence;
	std::chrono::steady_clock::time_point t_start, t_end;
	std::chrono::duration<double, std::milli> duration;
	
	start_time = continue_time + 1;  // start time can be 1 or any other number (less than duration) for restarts (repeats)
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
		// executes the simulation for each individual lattice site, allowing for parallelization
		SiteStepper stepper(*this);

		// index lattice sites
		int i, j;
		#pragma omp single
		{
			t_end = std::chrono::steady_clock::now();
			duration = (t_end - t_start) / 1000.;
			fprintf(stdout, "Done setting up RNGs in %.4f seconds\n", duration.count());
			fprintf(stdout, "Starting at time step %d\n", start_time);
		}
		// iterate through time steps (each thread goes through same time steps)
		for (int time_step = start_time; time_step < max_time_step + 1; time_step++) {
			// set start time for time step
			#pragma omp single
			t_start = std::chrono::steady_clock::now();

			// parallelized: distribute sites in lattice among threads
			// must be monotonic--site index must only ever increase because of discardRandom
			#pragma omp for schedule(monotonic:static)
			for (int index = 0; index < lattice_size * lattice_size; index++) {
				// lattice index to 2D index
				i = index / lattice_size;
				j = index - i * lattice_size;
				// update single site in 'next_lattice' and 'next_dispersal_lattice'
				stepper.updateSingleSite(i, j, time_step);
			}
			// one thread calculates persistence
			#pragma omp single
			persistence = getPersistence();  // get the number of species currently coexisting in the lattice
			// if persistence is below the minimum allowed persistence, the simulation is reinitialized with new random seeds
			if (persistence < min_persistence) {
				#pragma omp single
				{
					fprintf(stdout, "\n\nToo many extinctions, reinitializing\n\n");
					reinitializeSimulation(time_step);
					// save new initial state
					saveCompetition();
					saveLattice(0);
					persistence = 0;
				}
				// pass new seeds to each thread's SiteStepper RNG
				stepper.initializeRandomGenerator();
				time_step = 0;
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

int Ecolattice::getPersistence(void) {
	/* calculate number of species in lattice given only the species_abundance vector, 
		i.e., does not go through lattice. */
	int p = 0;
	for (int i = 0; i < num_species; i++) {
		if (species_abundance[i] != 0)
			p++;
		if (species_abundance[i] < 0) {
			fprintf(stderr, "Error, negative abundance for species %d\n", i + 1);
 			exit(-1);
		}
	}
	return p;
}

void Ecolattice::nextToThis(void) {
	/* this method updates the species in the lattice at this time step to the next time step,
	and sets the next time step to be unoccupied. does the same for seeds in the dispersal
	lattice */

	#pragma omp for
	for (int i = 0; i < lattice_size; i++) {
		for (int j = 0; j < lattice_size; j++) {
			lattice[i][j] = next_lattice[i][j];
			next_lattice[i][j] = 0;
			for (int k = 0; k < num_species; k++) {
				dispersal_lattice[i][j][k] = next_dispersal_lattice[i][j][k];
				next_dispersal_lattice[i][j][k] = 0;
			}
		}
	}
	return;
}


void Ecolattice::reinitializeSimulation(int time_step) {
	/* used to reinitialize the simulation after the number of species drops below the minimum persistence.
		method deletes previously written output files, draws new random seeds, and initializes simulation. */

	// delete previously written output files if simulation is being reinitialized after persistence is too low
	if (id == 0) {
		std::string fname;
			fname = outfile_base + "_competition.csv";
			remove(fname.c_str());
		for (int i = 0; i < time_step; i++) {
			fname = outfile_base+"_" + std::to_string(i) + ".csv";
			remove(fname.c_str());
		}
		for (int i = 1; i < num_species + 1; i++) { 
			fname = outfile_base + "_dispersal_s" + std::to_string(i) + "_0.csv";
			remove(fname.c_str());
			fname = outfile_base + "_dispersal_s" + std::to_string(i) + "_1.csv";
			remove(fname.c_str());
		}
		num_restarts++;
	}

	// draw new random seeds and reset random count
	drawRandomSeeds();
	random_count = 0;

	if (competition_filename.size() == 0 || continue_time != 0)
		// new simulation
		initializeRandomSimulation();
	else if (continue_time == 0)
		// replicate simulation
		initializeReplicateSimulation();

	// calculate properties of pairwise and community-level interactions
	getSpeciesAbundance();
	getImbalanceMean();
	getDiscreteTransitivity();
	getFecundityGrowthCorrelation();
	// RNG is set to jump ahead to the maximum number of random draws that could have already been used by the RNG when initializing the simulation
	unsigned long long max_random_count = 1000. * (4. * num_species * num_species + 5. * lattice_size * (unsigned long long) lattice_size);

	if (random_count > max_random_count) {
		if (id == 0) {
			fprintf(stderr, "Error, too many random numbers used to generate initial conditions.\n");
			fprintf(stderr, "Probable causes are the parameterization of TNormal distribution or severe competition correlation\n");
		}
		exit(0);
	}
	discardRandom(max_random_count - random_count);
	random_count = 0;
	return;
}
