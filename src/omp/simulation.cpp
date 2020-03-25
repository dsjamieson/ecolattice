
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2018			
	 *
	 *	The Simulation class constructor is defined here.
	 *	requires an imput file containing the simulation 
	 *	parameters. public methods for getting, setting and
	 *	adding to values of private arrays and variable are 
	 *	here, as well as the methods for saving output.
	 *
	 ***********************************************************/

#include "simulation.hpp"

Simulation::Simulation(std::string filename, int p_id) {
	
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
	setRandomSeeds();
	// checks for any obvious format issues with input parameter file 
	checkInputFormat();

	getParameter(num_threads, "NumThreads", 0);
	if (num_threads < 0) {
		if (id == 0)
			fprintf(stderr, "Error, NumThreads must be positive\n");
		exit(0);
	}
	omp_set_num_threads(num_threads);
	fprintf(stdout, "Ecolattice running in parallel (OMP) with %d threads\n", num_threads);
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
	// if a competition file is provided, simulation will run as a replicate, with the same random parameters but new matrix initiation and Monte Carlo dynamics
	// if simulation is a restart, there must be a competition file provided
	if (continue_time != -1) {
		getParameter(competition_filename, "CompetitionFile", 1);
	}
	else {
		getParameter(competition_filename, "CompetitionFile", 0);
	}

	// parameters read from file (argument: filename)
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

	if (competition_filename.size() != 0) {
		if(continue_time == 0) {
			// rerun simulation with pre-defined seeds (restart simulation)
			fprintf(stdout, "Initializing restarted simulation\n");
			loadSeeds();
			initializeRandomSimulation();
		}
		else if (continue_time > 0) {
			// continue a previous failed simulation or extend a previous simulation (continue simulation)
			fprintf(stdout, "Initialing continuation of previous simulation\n");
			initializeContinueSimulation();
		}
		else {
			// run simulation with pre-defined competition parameters but new initial conditions (replicate simulation)
			continue_time = 0;
			fprintf(stdout, "Initializing a replicate simulation\n");
			initializeReplicateSimulation();
		}
	}
	else {
		// start a new simulation
		continue_time = 0;
		fprintf(stdout, "Initializing new simulation\n");
		initializeRandomSimulation();
	}

	// calculate properties of pairwise and community-level interactions
	getSpeciesAbundance();
	getImbalanceMean();
	getDiscreteTransitivity();
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

	getFecundityGrowthCorrelation();

	// RNG is set to jump ahead to the maximum number of random draws that could have already been used by the RNG when initializing the simulation
	max_random_count = 1000. * (4. * num_species * num_species + 5. * lattice_size * (unsigned long long) lattice_size);

	if (random_count > max_random_count) {
		if (id == 0) {
			fprintf(stderr, "Error, too many random numbers used to generate initial conditions.\n");
			fprintf(stderr, "Probable causes are the parameterization of TNormal distribution or severe competition correlation\n");
		}
		exit(0);
	}
	discardRandom(max_random_count - random_count);
	random_count = 0;  // reset random count before simulation begins
	t_end = std::chrono::steady_clock::now();
	duration = (t_end - t_start) / 1000.;
	fprintf(stdout, "Done initialization in %.4f seconds\n", duration.count());
}


void Simulation::initializeRandomSimulation() {
	/* initializes the simulation lattice with species locations, and draws random variates for species-specific parameters
	(dispersal, competition, etc.). also checks that parameter values are appropriate. */

	// send random seeds to RNG
	seedGenerator();

	// set species specific parameters, potentially random
	getParameter(delta, "Delta", 2);
	for (int i = 0; i < num_species; i++) {
		neighborhood_size[i] = (2 * delta[i] + 1) * (2 * delta[i] + 1) - 1;
	}
	setRandomParameter(species_occupancy, "SpeciesOccupancy", 3);
	setRandomParameter(juvenile_survival_probability, "JuvenileSurvival", 2);
	setRandomParameter(adult_survival_probability, "AdultSurvival", 2);
	setRandomParameter(maximum_competition, "MaximumCompetition", 2);
	setRandomParameter(dispersal_probability, "DispersalProbability", 2);
	setRandomParameter(dispersal_length, "DispersalLength", 4);
	for (int i = 0; i < num_species; i++) {
		if ((int) round(2 * dispersal_length[i]) > lattice_size) {
			if (id == 0)
				fprintf(stderr, "Error, DispersalLength must be less than half of LatticeSize\n");
			exit(0);
		}
	}
	setRandomParameter(intrinsic_fecundity, "Fecundity", 4);
	
	// set competition parameters
	// competition lower and upper bound
	getParameter(competition_lower_bound, "CompetitionLower", 0);
	if (fabs(competition_lower_bound) > 1.) {
		if (id == 0)
			fprintf(stderr, "Error, CompetitionLower must be between -1 and 1\n");
		exit(0);
	}
	getParameter(competition_upper_bound, "CompetitionUpper", 0);
	if (fabs(competition_upper_bound) > 1.) {
		if (id == 0)
			fprintf(stderr, "Error, CompetitionUpper must be between -1 and 1\n");
		exit(0);
	}
	competition_mean = (competition_lower_bound + competition_upper_bound) / 2.;
	// competition for diagonal elements can have different parameters
	getParameter(competition_diag_lower_bound, "CompetitionDiagLower", 0);
	if (fabs(competition_diag_lower_bound) > 1.) {
		if (id == 0)
			fprintf(stderr, "Error, CompetitionDiagLower must be between -1 and 1\n");
		exit(0);
	}
	getParameter(competition_diag_upper_bound, "CompetitionDiagUpper", 0);
	if (fabs(competition_diag_upper_bound) > 1.) {
		if (id == 0)
			fprintf(stderr, "Error, CompetitionDiagUpper must be between -1 and 1\n");
		exit(0);
	}
	competition_diag_mean = (competition_diag_lower_bound + competition_diag_upper_bound) / 2.;

	// competition type (uniform, truncated normal), mean, and std. deviation
	getParameter(competition_type, "CompetitionType", 0);
	if (competition_type.compare("TNormal") == 0) {
		getParameter(competition_mean, "CompetitionMean", 0);
		if (competition_mean < competition_lower_bound || competition_mean > competition_upper_bound) {
			if (id == 0)
				fprintf(stderr, "Error, CompetitionMean must be between CompetitionLower and CompetitionUpper\n");
			exit(0);
		}
		getParameter(competition_diag_mean, "CompetitionDiagMean", 0);
		if (competition_diag_mean < competition_diag_lower_bound || competition_diag_mean > competition_diag_upper_bound) {
			if (id == 0)
				fprintf(stderr, "Error, CompetitionDiagMean must be between CompetitionDiagLower and CompetitionDiagUpper\n");
			exit(0);
		}
		getParameter(competition_sdev, "CompetitionSdev", 1);
		getParameter(competition_diag_sdev, "CompetitionDiagSdev", 1);
	}
	// optional competition structural features
	// correlation between growth and fecundity competition
	getParameter(competition_correlation, "CompetitionCorr", 0);
	if (fabs(competition_correlation) > 1.) {
		if (id == 0)
			fprintf(stderr, "Error, CompetitionCorr must be between -1 and 1\n");
		exit(0);
	}
	// imbalance between species pairs in competition matrices
	getParameter(imbalance, "Imbalance", 0);
	if (imbalance < 0 || imbalance > 1) {
		if (id == 0)
			fprintf(stderr, "Error, Imbalance must be between 0 and 1\n");
		exit(0);
	}
	// transitivity of matrices, 0 = random, 1 = max. transitivity, -1 = max intransitivity
	// relative hierarchy between fecundity/growth competition, +/-1 = equal/inverted, 0 = random 
	getParameter(fecundity_transitivity_type, "FecundityTransitivity", 0);
	if (fabs(fecundity_transitivity_type) != 1. && fecundity_transitivity_type != 0.) {
		if (id == 0)
			fprintf(stderr, "Error, current implementation only allows for FecundityTransitivity 0 (random), maximum (1), or minimum (-1)\n");
		exit(0);
	}
	getParameter(growth_transitivity_type, "GrowthTransitivity", 0);
	if (fabs(growth_transitivity_type) != 1. && growth_transitivity_type != 0.) {
		if (id == 0)
			fprintf(stderr, "Error, current implementation only allows for GrowthTransitivity 0 (random), maximum (1), or minimum (-1)\n");
		exit(0);
	}
	getParameter(fecundity_growth_relative_hierarchy, "RelativeHierarchy", 0);
	if (fabs(fecundity_growth_relative_hierarchy) != 1. && fecundity_growth_relative_hierarchy != 0.) {
		if (id == 0)
			fprintf(stderr, "Error, RelativeHierarchy must be +/-1 (equal/inverted), or 0 (independent)\n");
		exit(0);
	}
	if (fecundity_growth_relative_hierarchy != 0 && (fecundity_transitivity_type == 0 || growth_transitivity_type == 0)) {
		if (id == 0)
			fprintf(stderr, "Error, if RelativeHierarchy is not zero, neither FecundityTransitivity nor GrowthTransitivity can be zero\n");
		exit(0);
	}
	if (fecundity_growth_relative_hierarchy != 0 && (fecundity_transitivity_type != growth_transitivity_type)) {
		if (id == 0)
			fprintf(stderr, "Error, if RelativeHierarchy is not zero, FecundityTransitivity and GrowthTransitivity must be equal\n");
		exit(0);
	}
	if (fecundity_growth_relative_hierarchy == -1. && (fecundity_transitivity_type == -1. && growth_transitivity_type == -1.)) {
		if (id == 0)
			fprintf(stderr, "Error, if both FecundityTransitivity and GrowthTransitivity are -1 (intransitive), RelativeHierarchy cannot be -1 (inverted).\n");
		exit(0);
	}
	if ((competition_correlation != 0) + (imbalance != 0.5) + ((fabs(fecundity_transitivity_type) + fabs(growth_transitivity_type)) != 0) > 1) {
		if (id == 0)
			fprintf(stderr, "Error, only one of CompetitionCorr, Imbalance, and (Fecundity/Growth)Transitivity can be set\n");
		exit(0);
	}

	// initializate competition
	// first, draw random variates from uniform or truncated normal distribution to fill competition matrices
	if (competition_type.compare("Uniform") == 0) {
		// if growth and fecundity matrices are correlated or anticorrelated (growth/reproduction trade-off)
		if (competition_correlation != 0) 
			initializeUniformCorrelatedCompetition();
		else
			initializeUniformCompetition();
	}
	else if (competition_type.compare("TNormal") == 0){
		// if growth and fecundity matrices are correlated or anticorrelated (growth/reproduction trade-off)
		if (competition_correlation != 0)
			initializeTNormalCorrelatedCompetition();
		else
			initializeTNormalCompetition();
	}
	else {
		if (id == 0)
			fprintf(stderr, "CompetitionType must be either Uniform or TNormal\n");
		exit(0);
	}
	// if competition is imbalanced (one species has a stronger effect on the second than the second has on the first
	if (imbalance != 0.5)
		imbalanceCompetition();
	// if competition is transitive (community forms a hierarchy) or intransitive (community has loops breaking hierarchy)
	if (fecundity_transitivity_type != 0 || growth_transitivity_type!=0)
		setCompetitionTransitivity();

	// RNG discard after drawing random parameters, as some use rejection sampling
	unsigned long long max_random_count = 1000. * (4. * num_species * (unsigned long long) num_species);

	if (random_count > max_random_count) {
		if (id == 0) {
			fprintf(stderr, "Error, too many random numbers used to generate competition and parameters.\n");
			fprintf(stderr, "Probable causes are the parameterization of TNormal distribution or severe competition correlation\n");
		}
		exit(0);
	}
	discardRandom(max_random_count - random_count);

	initializeLattice();

	return;
}

void Simulation::initializeReplicateSimulation() {
	/* method used to start a new simulation with the same random parameters, used for replicates.
	uses competition matrices, fecundities, occupancies, etc. from file, specified in 
	"CompetitionFile." replicates have the same parameters, but are different realizations
	(i.e., dynamics will differ). */

	// send seeds to RNG
	seedGenerator();
	// read in Delta from parameter file and all other parameters from competition file
	getParameter(delta, "Delta", 2);
	for (int i = 0; i < num_species; i++) {
		neighborhood_size[i] = (2 * delta[i] + 1) * (2 * delta[i] + 1) - 1;
	}
	loadCompetition();
	// RNG discard after loading random parameters, as some used rejection sampling
	unsigned long long max_random_count = 1000. * (4. * num_species * (unsigned long long) num_species);

	if (random_count > max_random_count) {
		if (id == 0) {
			fprintf(stderr, "Error, too many random numbers used to generate competition and parameters.\n");
			fprintf(stderr, "Probable causes are the parameterization of TNormal distribution or severe competition correlation\n");
		}
		exit(0);
	}
	discardRandom(max_random_count - random_count);
	initializeLattice();
	return;
}

void Simulation::initializeContinueSimulation() {
	/* method used if a previous simulation failed before completing or if you want to extend the simulation.
	continues the simulation from where it left off. this method initializes the lattice, reloads the parameters from
	the previous simulation, and starts up the RNG for the appropriate time step. restarted simulations are identical
	in both parameters and dynamics. */

	// load seeds from competition file and sends to RNG
	loadSeeds();
	seedGenerator();
	// load dispersal from file
	loadDispersal();
	// read in Delta from parameter file and all other parameters from competition file
	getParameter(delta, "Delta", 2);
	for (int i = 0; i < num_species; i++) {
		neighborhood_size[i] = (2 * delta[i] + 1) * (2 * delta[i] + 1) - 1;
	}
	loadCompetition();
	loadLattice();
	return;
}

void Simulation::reinitializeSimulation(int time_step) {
	/* method used for two purposes: 1) to reinitialize worker copies of the simulation given the random seeds drawn
	by the central task, and 2) to reinitialize the simulation after the number of species drops below the minimum persistence.
	method deletes already written output files (only for purpose 2) */

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
	random_count = 0;

	if (competition_filename.size() == 0 || continue_time != 0)
		// new simulation, seeds from id = 0
		initializeRandomSimulation();
	else if (continue_time == 0)
		// replicate simulation, seeds from id = 0
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

void Simulation::allocSimulation() {
	/* allocate memory for all arrays used in simulations, including parameter arrays, lattice with species locations, and dispersal lattice with seed locations. */

	species_abundance.resize(num_species, 0);
	delta.resize(num_species, 0);
	neighborhood_size.resize(num_species, 0);
	juvenile_survival_probability.resize(num_species, 0);
	adult_survival_probability.resize(num_species, 0);
	maximum_competition.resize(num_species, 0);
	dispersal_probability.resize(num_species, 0);
	fecundity_rank.resize(num_species, 0);
	growth_rank.resize(num_species, 0);
	species_occupancy.resize(num_species, 0);
	dispersal_length.resize(num_species, 0);
	intrinsic_fecundity.resize(num_species, 0);
	competition_fecundity.resize(num_species);
	competition_growth.resize(num_species);
	fecundity_transitivity.resize(num_species);
	growth_transitivity.resize(num_species);
	for (int i = 0; i < num_species ; i++) {
		competition_fecundity[i].resize(num_species, 0);
		competition_growth[i].resize(num_species, 0);
		fecundity_transitivity[i].resize(num_species, 0);
		growth_transitivity[i].resize(num_species, 0);
	}

	lattice.resize(lattice_size);
	next_lattice.resize(lattice_size);
	dispersal_lattice.resize(lattice_size);
	next_dispersal_lattice.resize(lattice_size);
	for (int i = 0; i < lattice_size; i++) {
		lattice[i].resize(lattice_size, 0);
		next_lattice[i].resize(lattice_size, 0);
		dispersal_lattice[i].resize(lattice_size);
		next_dispersal_lattice[i].resize(lattice_size);
		for (int j = 0; j < lattice_size; j++) {
			lattice[i][j] = 0;
			next_lattice[i][j] = 0;
			dispersal_lattice[i][j].resize(num_species, 0);
			next_dispersal_lattice[i][j].resize(num_species, 0);
		}
	}
	return;
}


void Simulation::initializeLattice() {
	/* initializes lattice with randomly located species, depending on the occupancy probability, which defines both the total
	occupancy of the lattice and the species specific probabilities. also determines whether individuals are juveniles or adults
	with equal probability. */

	std::bernoulli_distribution stage_dist(0.5);
	std::bernoulli_distribution occupy_dist(initial_occupancy);
	std::discrete_distribution<int> species_dist(species_occupancy.begin(), species_occupancy.end());

	for (int i = 0; i < lattice_size; i++) {
		for (int j = 0; j < lattice_size; j++) {
			lattice[i][j] = occupy_dist(generateRandom());
			if (stage_dist(generateRandom()) == 0 && lattice[i][j] != 0)
				lattice[i][j] *= -1;
			if (lattice[i][j] != 0) {
				lattice[i][j] *= 1 + species_dist(generateRandom());
			}
		}
	}
	return;
}

void Simulation::getSpeciesAbundance() {
	/* calculate the number of individuals of each species in the lattice. This only happens after initialization, before
 * 	   any simulation time steps. Subsequently, species abundances are kept track of using the incrementSpeciesAbundance
 * 	   and decrementSpeciesAbundance methods. */

	for (int i = 0; i < num_species; i++)
		species_abundance[i] = 0;

	#pragma omp parallel for
	for (int i = 0; i < lattice_size; i++) {
		for (int j = 0; j < lattice_size; j++) {
			int s = abs(lattice[i][j]);
			if (s != 0) {
				#pragma omp atomic
				species_abundance[s - 1]++;
			}
		}
	}
	return;
}


void Simulation::resetLattice() {
	/* for the current time step, set all elements of the lattice and the dispersal lattice to 0. used in simulation for the first time step (t = 0). */

	for (int i = 0; i < lattice_size; i++) {
		for (int j = 0; j < lattice_size; j++) {
			lattice[i][j] = 0;
			for (int k = 0; k < num_species; k++) {
				dispersal_lattice[i][j][k] = 0.;
			}
		}
	}
	return;
}


void Simulation::resetNextLattice() {
	/* for the next time step, set all elements of the lattice and the dispersal lattice to 0. used in simulation so that workers
	can reset their local copies of 'next_lattice' and 'next_dispersal_lattice' */

	for (int i = 0; i < lattice_size; i++) {
		for (int j = 0; j < lattice_size; j++) {
			next_lattice[i][j] = 0;
			for (int k = 0; k < num_species; k++) {
				next_dispersal_lattice[i][j][k] = 0.;
			}
		}
	}
	return;
}


std::mt19937& Simulation::generateRandom() {
	/* add to the random count and get a random draw from the global RNG. */

	random_count += 2;
	return global_random_generator;
}


void Simulation::discardRandom(unsigned long long n) {
	/* add to the random count and discard values from the global RNG. */

	global_random_generator.discard(n);
	random_count += n;
	return;
}

void Simulation::nextToThis() {
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


void Simulation::saveLattice(int time_step) {
	/* saves the species locations in the lattice from the current time step to file */

	std::ofstream lattice_file;
	lattice_file.open(outfile_base+"_" + std::to_string(time_step) + ".csv", std::ios::out | std::ios::trunc);

	if (!lattice_file.is_open()) {
			if (id == 0)
				fprintf(stderr, "Error, could not open time step %d lattice file file to load\n", time_step);
			exit(0);
	}
	for (int i = 0; i < lattice_size; i++) {
		for (int j = 0; j < lattice_size; j++) {
			if (lattice[i][j] < 0)
				lattice_file << lattice[i][j];
			else
				lattice_file << " " << lattice[i][j];
			if (j != lattice_size - 1)
				lattice_file << ", ";
			}
			lattice_file << std::endl;
		}

	lattice_file.close();

	return;
}


void Simulation::saveCompetition() {
	/* before the first time step, saves the parameters from this simulation to file, including species
	occupancy, survival, maximum competition, dispersal probability, dispersal length, intrinsic fecundity
	fecundity and growth competition matrices. also saves metrics calculated from the competition matrices,
	i.e., transitivity. */

	int i, j;

	std::ofstream competition_file;
	competition_file.open(outfile_base + "_competition.csv", std::ios::out | std::ios::trunc);

	if (!competition_file.is_open()) {
			if (id == 0)
				fprintf(stderr, "Error, could not open competition file to save\n");
			exit(0);
	}

	competition_file << "# Seeds:" << std::endl;
	for (i = 0; i < 5; i++) {
		competition_file << " " << seeds[i];
		if (i < 4)
			competition_file << ",";
	} 
	competition_file << std::endl;

	competition_file << "# Species Occupancy:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << species_occupancy[i];
		if (i < num_species-1)
			competition_file << ",";
	}
	competition_file << std::endl;

	competition_file << "# Juvenile Survival:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << juvenile_survival_probability[i];
		if (i < num_species-1)
			competition_file << ",";

	}
	competition_file << std::endl;

	competition_file << "# Adult Survival:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << adult_survival_probability[i];
		if (i < num_species-1)
			competition_file << ",";

	}
	competition_file << std::endl;

	competition_file << "# Maximum Competition:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << maximum_competition[i];
		if (i < num_species-1)
			competition_file << ",";
	}
	competition_file << std::endl;

	competition_file << "# Dispersal probability:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << dispersal_probability[i];
		if (i < num_species-1)
			competition_file << ",";
	}
	competition_file << std::endl;

	competition_file << "# Dispersal length:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << dispersal_length[i];
		if (i < num_species-1)
			competition_file << ",";
	}
	competition_file << std::endl;

	competition_file << "# Intrinsic fecundity:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << intrinsic_fecundity[i];
		if (i < num_species-1)
			competition_file << ",";
	}
	competition_file << std::endl;
	
	competition_file << "# Fecundity competition:" << std::endl;
	for (i = 0; i < num_species; i++) {
		for (j = 0; j < num_species; j++) {
			if (competition_fecundity[i][j] < 0)
				competition_file << competition_fecundity[i][j];
			else
				competition_file << " " << competition_fecundity[i][j];
			if (j != num_species - 1)
				competition_file << ", ";
			}
			competition_file << std::endl;
		}

	competition_file << "# Growth competition:" << std::endl;
	for (i = 0; i < num_species; i++) {
		for (j = 0; j < num_species; j++) {
			if(competition_growth[i][j] < 0)
				competition_file << competition_growth[i][j];
			else
				competition_file << " " << competition_growth[i][j];
			if (j != num_species - 1)
				competition_file << ", ";
			}
			competition_file << std::endl;
		}

	competition_file << "# Fecundity transitivity ranks:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << fecundity_rank[i];
		if (i < num_species-1)
			competition_file << ",";
	}
	competition_file << std::endl;

	competition_file << "# Growth transitivity ranks:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << growth_rank[i];
		if (i < num_species-1)
			competition_file << ",";
	}
	competition_file << std::endl;

	competition_file << "# Fecundity transitivity:" << std::endl;
	for (i = 0; i < num_species; i++) {
		for (j = 0; j < num_species; j++) {
			if (fecundity_transitivity[i][j] < 0)
				competition_file << fecundity_transitivity[i][j];
			else
				competition_file << " " << fecundity_transitivity[i][j];
			if (j != num_species - 1)
				competition_file << ", ";
			}
			competition_file << std::endl;
		}

	competition_file << "# Growth transitivity:" << std::endl;
	for (i = 0; i < num_species; i++) {
		for (j = 0; j < num_species; j++) {
			if (growth_transitivity[i][j] < 0)
				competition_file << growth_transitivity[i][j];
			else
				competition_file << " " << growth_transitivity[i][j];
			if (j != num_species - 1)
				competition_file << ", ";
			}
			competition_file << std::endl;
		}

	competition_file <<  "# Fecundity imbalance mean:" << std::endl  << fecundity_imbalance_mean <<  std::endl;
	competition_file <<  "# Growth imbalance mean:" << std::endl  << growth_imbalance_mean <<  std::endl;
	competition_file <<  "# Fecundity relative intransitivity:" << std::endl  << fecundity_relative_intransitivity <<  std::endl;
	competition_file <<  "# Growth relative intransitivity:" << std::endl  << growth_relative_intransitivity <<  std::endl;
	competition_file <<  "# Fecundity-growth cross correlation:" << std::endl  << fecundity_growth_correlation <<  std::endl;
	competition_file <<  "# Number of restarts:" << std::endl  << num_restarts <<  std::endl;
	
	competition_file.close();

	return;

}

/* methods used to "get" and "set" values inside the main function. "get" refers to returning values
from the simulation object. "set" refers to changing values within the simulation given values in main. */

int Simulation::getLatticeSize() {
	return lattice_size;
}


int Simulation::getNumSpecies() {
	return num_species;
}

double Simulation::getGerminationProbability() {
	return germination_probability;
}

double Simulation::getAdultSurvivalProbability(int s) {
	return adult_survival_probability[s];
}

double Simulation::getJuvenileSurvivalProbability(int s) {
	return juvenile_survival_probability[s];
}

double Simulation::getCompetitionGrowth(int s1, int s2) {
	return competition_growth[s1][s2];
}

double Simulation::getCompetitionFecundity(int s1, int s2) {
	return competition_fecundity[s1][s2];
}

double Simulation::getDelta(int s) {
	return delta[s];
}

double Simulation::getDispersalLength(int s) {
	return dispersal_length[s];
}

double Simulation::getDispersalProbability(int s) {
	return dispersal_probability[s];
}

double Simulation::getIntrinsicFecundity(int s) {
	return intrinsic_fecundity[s];
}

double Simulation::getMaximumCompetition(int s) {
	return maximum_competition[s];
}

int Simulation::getMaxTimeStep() {
	return max_time_step;
}


int Simulation::getSite(int i, int j) {
	return lattice[i][j];
}


int Simulation::getNextSite(int i, int j) {
	return next_lattice[i][j];
}


unsigned int Simulation::getRandom() {
	random_count++;
	return global_random_generator();

}

unsigned long long Simulation::getRandomCount() {
	return random_count;
}


double Simulation::getDispersal(int i, int j, int s) {
	return dispersal_lattice[i][j][s];
}


double Simulation::getNextDispersal(int i, int j, int s) {
	return next_dispersal_lattice[i][j][s];
}


void Simulation::setSite(int i, int j, int s) {
	lattice[i][j] = s;
	return;
}

void Simulation::setNextSite(int i, int j, int s) {
	next_lattice[i][j] = s;
	return;
}

void Simulation::addSite(int i, int j, int s) {
	lattice[i][j] += s;
	return;
}

void Simulation::setDispersal(int i, int j, int s, double t) {
	dispersal_lattice[i][j][s] = t;
	return;
}

void Simulation::addDispersal(int i, int j, int s, double t) {
	dispersal_lattice[i][j][s] += t;
	return;
}


int Simulation::getPersistence() {

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

void Simulation::printSpeciesAbundance(void) {

	for (int i = 0; i < num_species; i++) {
		fprintf(stdout, "s%d: %d\n", i + 1, species_abundance[i]);
	}
	return;
}

void Simulation::printNextLattice(void) {

	for (int i = 0; i < lattice_size; i++) {
		for (int j = 0; j < lattice_size; j++)
			fprintf(stdout, "%d ", next_lattice[i][j]);
		fprintf(stdout, "\n");
	}
	return;
}


int Simulation::getMinPersistence() {
	return min_persistence;
}

unsigned int Simulation::getSeed(int i) {
	return seeds[i];
}

void Simulation::setSeed(int i, unsigned int s) {
	seeds[i] = s;
	return;
}
