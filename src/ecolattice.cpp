
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2020		
	 *
	 *	the Ecolattice class constructor is defined here.
	 *	requires an input file containing the simulation 
	 *	parameters.
	 *
	 ***********************************************************/

#include "ecolattice.hpp"

Ecolattice::Ecolattice(std::string filename, int p_id) {
	
	id = p_id;
	random_count = 0;

	// default parameters
	num_threads = 1;
	continue_time = -1;
	competition_type = "Uniform";
	competition_upper_bound = 0;
	competition_lower_bound = -1;
	competition_diag_upper_bound = 0;
	competition_diag_lower_bound = -1;
	competition_correlation = 0;
	imbalance = 0.5;
	fecundity_transitivity_type = 0.;
	growth_transitivity_type = 0.;
	fecundity_growth_relative_hierarchy = 0.;
	min_persistence = 0;
	num_restarts = 0;
	parameter_filename = filename;

	std::chrono::steady_clock::time_point t_start, t_end;
	std::chrono::duration<double, std::milli> duration;
	t_start = std::chrono::steady_clock::now();

	// draw random seeds used in simulation
	drawRandomSeeds();
	// checks for any obvious format issues with input parameter file 
	checkInputFormat();

	// parallelized OMP code can run on any number of threads
	getParameter(num_threads, "NumThreads", 0);
	if (num_threads < 1) {
		if (id == 0)
			fprintf(stderr, "Error, NumThreads must be positive\n");
		exit(0);
	}
#ifdef OMP
	omp_set_num_threads(num_threads);
	fprintf(stdout, "Ecolattice running in parallel (OMP) with %d threads\n", num_threads);
#else
	fprintf(stdout, "Ecolattice running in serial\n");
#endif
	/*#pragma omp parallel
	{
		int threads = omp_get_max_threads();
		#pragma omp single
		fprintf(stdout, "max_threads: %d\n", threads);
		exit(0);
	}*/

	// checks if this simulation is a continuation of a previous simulation given time step (parameter: continue_time)
	int temp_continue_time;
	if(getParameter(temp_continue_time, "ContinueTime", 0)) {
		if (temp_continue_time < 0) {
			if (id == 0)
				fprintf(stderr, "Error, ContinueTime must be non-negative\n");
			exit(0);
		}
		continue_time = temp_continue_time;
	}
	// checking if competition file is provided
	// will run either continue, replicate, or restarted simulation	
	if (continue_time != -1) {
		getParameter(competition_filename, "CompetitionFile", 1);
	}
	else {
		getParameter(competition_filename, "CompetitionFile", 0);
	}

	// checking parameters read from file (argument: filename)
	getParameter(lattice_size, "LatticeSize", 1);
	if (lattice_size < 0) {
		if (id == 0)
			fprintf(stderr, "Error, LatticeSize must be positive\n");
		exit(0);
	}
	getParameter(num_species, "Species", 1);
	getParameter(max_time_step, "MaxTimeStep", 1);
	if(continue_time >= max_time_step) {
		if (id == 0)
			fprintf(stderr, "Error, ContinueTime must be less than MaxTimeStep\n");
		exit(0);
	}
	getParameter(initial_occupancy, "InitialOccupancy", 1);
	if (initial_occupancy > 1 || initial_occupancy < 0) {
		if (id == 0)
			fprintf(stderr, "Error, InitialOccupancy must be between 0 and 1\n");
		exit(0);
	}
	getParameter(germination_probability, "GerminationProbability", 1);
	if (germination_probability > 1 || germination_probability < 0) {
		if (id == 0)
			fprintf(stderr, "Error, GerminationProbability must be between 0 and 1\n");
		exit(0);
	}
	getParameter(outfile_base, "OutfileBase", 1);
	getParameter(outfile_dir, "OutfileDir", 0);
	if (outfile_dir.size() != 0) {
		struct stat buf;
		if (stat(outfile_dir.c_str(), &buf) == 0)
			outfile_base = outfile_dir + "/" + outfile_base;
		else {
			if (id == 0)
				fprintf(stderr, "Error, no directory %s\n", outfile_dir.c_str());
			exit(0);
		}
	}
	getParameter(min_persistence, "MinPersistence", 0);
	if (min_persistence > num_species || min_persistence < 0) {
		if (id == 0)
			fprintf(stderr, "Error, MinPersistence must be between 0 and Species\n");
		exit(0);
	}

	// allocate simulation arrays
	allocSimulation();

	// different types of simulations can be initialized: restart, continue, replicate, or random (default)
	if (competition_filename.size() != 0) {
		if (continue_time == 0) {
			// rerun identical simulation with pre-defined seeds (restart simulation)
			// uses seeds from competition file
			fprintf(stdout, "Initializing restarted simulation\n");
			initialization_scheme = 1;
			loadSeeds();
			initializeRandomSimulation();
		}
		else if (continue_time > 0) {
			// continue a previous failed simulation or extend a previous simulation (continue simulation)
			// uses seeds from competition file
			fprintf(stdout, "Initializing continuation of previous simulation\n");
			initialization_scheme = 2;
			initializeContinueSimulation();
		}
		else {
			// run simulation with pre-defined competition parameters but new initial conditions (replicate simulation)
			continue_time = 0;
			fprintf(stdout, "Initializing a replicate simulation\n");
			initialization_scheme = 3;
			initializeReplicateSimulation();
		}
	}
	else {
		// start a new simulation
		continue_time = 0;
		fprintf(stdout, "Initializing new simulation\n");
		initialization_scheme = 0;
		initializeRandomSimulation();
	}

	// calculate properties of pairwise and community-level interactions
	getSpeciesAbundance();
	getImbalanceMean();
	getDiscreteTransitivity();
	getFecundityGrowthCorrelation();

	// checks for failure, e.g., does calculated transitivity match specified transitivity
	if (fecundity_relative_intransitivity != 1. && fecundity_transitivity_type == -1.) {
		fprintf(stderr, "Error, setting transitivity type failed\n");
		exit(0);
	}
	else if (fecundity_relative_intransitivity != 0. && fecundity_transitivity_type == 1.) {
		fprintf(stderr, "Error, setting transitivity type failed\n");
		exit(0);
	}
	if (growth_relative_intransitivity != 1. && growth_transitivity_type == -1.) {
		fprintf(stderr, "Error, setting transitivity type failed\n");
		exit(0);
	}
	else if (growth_relative_intransitivity != 0. && growth_transitivity_type == 1.) {
		fprintf(stderr, "Error, setting transitivity type failed\n");
		exit(0);
	}

	// RNG jumps ahead to the maximum number of random draws that could have already been used by the RNG when initializing
	max_random_count = 1000. * (4. * num_species * num_species + 5. * lattice_size * (unsigned long long) lattice_size);
	if (random_count > max_random_count) {
		if (id == 0) {
			fprintf(stderr, "Error, too many random numbers used to generate initial conditions.\n");
			fprintf(stderr, "Probable causes are the parameterization of TNormal distribution or severe competition correlation\n");
		}
		exit(0);
	}
	discardRandom(max_random_count - random_count);
	// reset random count before simulation begins
	random_count = 0;
	t_end = std::chrono::steady_clock::now();
	duration = (t_end - t_start) / 1000.;
	fprintf(stdout, "Done initialization in %.4f seconds\n", duration.count());
}
