
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

#include "simulation.h"

Simulation::Simulation(std::string filename, int p_id) {
	
	id = p_id;
	random_count = 0;

	// default parameters
	restart_time = 0;
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

	// draw random seeds used in simulation
	setRandomSeeds();
	// checks for any obvious format issues with input parameter file 
	checkInputFormat();

	// checks if this simulation is a continuation of a previous simulation given time step (parameter: restart_time)
	getParameter(&restart_time, "RestartTime", 0);
	if (restart_time < 0) {
		if (id == 0)
			fprintf(stderr, "Error, RestartTime must be positive\n");
		MPI_Finalize();
		exit(0);
	}

	// if a competition file is provided, simulation will run as a replicate, with the same random parameters but new matrix initiation and Monte Carlo dynamics
	// if simulation is a restart, there must be a competition file provided
	if (restart_time != 0) {
		getParameter(&competition_filename, "CompetitionFile", 1);
	}
	else {
		getParameter(&competition_filename, "CompetitionFile", 0);
	}

	// parameters read from file (argument: filename)
	getParameter(&lattice_size, "LatticeSize", 1);
	if (lattice_size < 0) {
		if (id == 0)
			fprintf(stderr, "Error, LatticeSize must be positive\n");
		MPI_Finalize();
		exit(0);
	}
	getParameter(&num_species, "Species", 1);
	getParameter(&max_time_step, "MaxTimeStep", 1);
	getParameter(&initial_occupancy, "InitialOccupancy", 1);
	if (initial_occupancy > 1 || initial_occupancy < 0) {
		if (id == 0)
			fprintf(stderr, "Error, InitialOccupancy must be between 0 and 1\n");
		MPI_Finalize();
		exit(0);
	}
	getParameter(&germination_probability, "GerminationProbability", 1);
	if (germination_probability > 1 || germination_probability < 0) {
		if (id == 0)
			fprintf(stderr, "Error, GerminationProbability must be between 0 and 1\n");
		MPI_Finalize();
		exit(0);
	}
	getParameter(&outfile_base, "OutfileBase", 1);
	getParameter(&outfile_dir, "OutfileDir", 0);
	if (outfile_dir.size() != 0) {
		struct stat buf;
		if (stat(outfile_dir.c_str(), &buf) == 0)
			outfile_base = outfile_dir + "/" + outfile_base;
		else {
			if (id == 0)
				fprintf(stderr, "Error, no directory %s\n", outfile_dir.c_str());
			MPI_Finalize();
			exit(0);
		}
	}
	getParameter(&min_persistence, "MinPersistence", 0);
	if (min_persistence > num_species || min_persistence < 0) {
		if (id == 0)
			fprintf(stderr, "Error, MinPersistence must be between 0 and Species\n");
		MPI_Finalize();
		exit(0);
	}

	// allocate simulation arrays
	allocSimulation();

	if (restart_time != 0) {
		// continue a previous failed simulation or extend a previous simulation (repeat simulation)
		initializeRestartSimulation();
	}
	else if (competition_filename.size() != 0) {
		// run simulation with pre-defined competition parameters but new initial conditions (replicate simulation)
		initializeReplicateSimulation();
	}
	else {
		// start a new simulation
		initializeRandomSimulation();
	}

	// calculate properties of pairwise and community-level interactions
	getSpeciesAbundance();
	getImbalanceMean();
	getDiscreteTransitivity();
	getFecundityGrowthCorrelation();

	// RNG is set to jump ahead to the maximum number of random draws that could have already been used by the RNG when initializing the simulation
	unsigned long long max_random_count =  (unsigned long long) 1000. * (4. * num_species * num_species + 5. * lattice_size * lattice_size);
	if (random_count > max_random_count) {
		if (id == 0)
			fprintf(stderr, "Error, too many random numbers used to generate initial conditions.\n");
			fprintf(stderr, "Probable causes are the parameterization of TNormal distribution or severe competition correlation\n");
		MPI_Finalize();
		exit(0);
	}
	discardRandom(max_random_count - random_count);
	random_count = 0;  // reset random count before simulation begins
}


void Simulation::initializeRandomSimulation() {
	/* initializes the simulation lattice with species locations, and draws random variates for species-specific parameters
	(dispersal, competition, etc.). also checks that parameter values are appropriate. */

	// send random seeds to RNG
	seedGenerator();

	// set species specific parameters, potentially random
	getParameter(delta, num_species, "Delta", 2);
	setRandomParameter(species_occupancy, num_species, "SpeciesOccupancy", 3);
	setRandomParameter(juvenile_survival_probability, num_species, "JuvenileSurvival", 2);
	setRandomParameter(adult_survival_probability, num_species, "AdultSurvival", 2);
	setRandomParameter(maximum_competition, num_species, "MaximumCompetition", 2);
	setRandomParameter(dispersal_probability,  num_species,"DispersalProbability", 2);
	setRandomParameter(dispersal_length, num_species, "DispersalLength", 4);
	setRandomParameter(intrinsic_fecundity, num_species, "Fecundity", 4);
	
	// set competition parameters
	// competition lower and upper bound
	getParameter(&competition_lower_bound, "CompetitionLower", 0);
	if (fabs(competition_lower_bound) > 1.) {
		if (id == 0)
			fprintf(stderr, "Error, CompetitionLower must be between -1 and 1\n");
		MPI_Finalize();
		exit(0);
	}
	getParameter(&competition_upper_bound, "CompetitionUpper", 0);
	if (fabs(competition_upper_bound) > 1.) {
		if (id == 0)
			fprintf(stderr, "Error, CompetitionUpper must be between -1 and 1\n");
		MPI_Finalize();
		exit(0);
	}
	competition_mean = (competition_lower_bound + competition_upper_bound) / 2.;
	// competition for diagonal elements can have different parameters
	getParameter(&competition_diag_lower_bound, "CompetitionDiagLower", 0);
	if (fabs(competition_diag_lower_bound) > 1.) {
		if (id == 0)
			fprintf(stderr, "Error, CompetitionDiagLower must be between -1 and 1\n");
		MPI_Finalize();
		exit(0);
	}
	getParameter(&competition_diag_upper_bound, "CompetitionDiagUpper", 0);
	if (fabs(competition_diag_upper_bound) > 1.) {
		if (id == 0)
			fprintf(stderr, "Error, CompetitionDiagUpper must be between -1 and 1\n");
		MPI_Finalize();
		exit(0);
	}
	competition_diag_mean = (competition_diag_lower_bound + competition_diag_upper_bound) / 2.;

	// competition type (uniform, truncated normal), mean, and std. deviation
	getParameter(&competition_type, "CompetitionType", 0);
	if (competition_type.compare("TNormal") == 0) {
		getParameter(&competition_mean, "CompetitionMean", 0);
		if (competition_mean < competition_lower_bound || competition_mean > competition_upper_bound) {
			if (id == 0)
				fprintf(stderr, "Error, CompetitionMean must be between CompetitionLower and CompetitionUpper\n");
			MPI_Finalize();
			exit(0);
		}
		getParameter(&competition_diag_mean, "CompetitionDiagMean", 0);
		if (competition_diag_mean < competition_diag_lower_bound || competition_diag_mean > competition_diag_upper_bound) {
			if (id == 0)
				fprintf(stderr, "Error, CompetitionDiagMean must be between CompetitionDiagLower and CompetitionDiagUpper\n");
			MPI_Finalize();
			exit(0);
		}
		getParameter(&competition_sdev, "CompetitionSdev", 1);
		getParameter(&competition_diag_sdev, "CompetitionDiagSdev", 1);
	}
	// optional competition structural features
	// correlation between growth and fecundity competition
	getParameter(&competition_correlation, "CompetitionCorr", 0);
	if (fabs(competition_correlation) > 1.) {
		if (id == 0)
			fprintf(stderr, "Error, CompetitionCorr must be between -1 and 1\n");
		MPI_Finalize();
		exit(0);
	}
	// imbalance between species pairs in competition matrices
	getParameter(&imbalance, "Imbalance", 0);
	if (imbalance < 0 || imbalance > 1) {
		if (id == 0)
			fprintf(stderr, "Error, Imbalance must be between 0 and 1\n");
		MPI_Finalize();
		exit(0);
	}
	// transitivity of matrices, 0 = random, 1 = max. transitivity, -1 = max intransitivity
	// relative hierarchy between fecundity/growth competition, +/-1 = equal/inverted, 0 = random 
	getParameter(&fecundity_transitivity_type, "FecundityTransitivity", 0);
	if (fabs(fecundity_transitivity_type) != 1. && fecundity_transitivity_type != 0.) {
		if (id == 0)
			fprintf(stderr, "Error, current implementation only allows for FecundityTransitivity 0 (random), maximum (1), or minimum (-1)\n");
		MPI_Finalize();
		exit(0);
	}
	getParameter(&growth_transitivity_type, "GrowthTransitivity", 0);
	if (fabs(growth_transitivity_type) != 1. && growth_transitivity_type != 0.) {
		if (id == 0)
			fprintf(stderr, "Error, current implementation only allows for GrowthTransitivity 0 (random), maximum (1), or minimum (-1)\n");
		MPI_Finalize();
		exit(0);
	}
	getParameter(&fecundity_growth_relative_hierarchy, "RelativeHierarchy", 0);
	if (fabs(fecundity_growth_relative_hierarchy) != 1. && fecundity_growth_relative_hierarchy != 0.) {
		if (id == 0)
			fprintf(stderr, "Error, RelativeHierarchy must be +/-1 (equal/inverted), or 0 (unrelated)\n");
		MPI_Finalize();
		exit(0);
	}
	if (fecundity_growth_relative_hierarchy != 0 && (fecundity_transitivity_type == 0 || growth_transitivity_type == 0)) {
		if (id == 0)
			fprintf(stderr, "Error, if RelativeHierarchy is not zero, neither FecundityTransitive nor GrowthTransitivity can be zero\n");
		MPI_Finalize();
		exit(0);
	}
	if ((competition_correlation!=0) + (imbalance!=0.5) + ((fabs(fecundity_transitivity_type) + fabs(growth_transitivity_type)) != 0) > 1) {
		if (id == 0)
			fprintf(stderr, "Error, only one of CompetitionCorr, Imbalance, and (Fecundity/Growth)Transitivity can be set\n");
		MPI_Finalize();
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
		MPI_Finalize();
		exit(0);
	}
	// if competition is imbalanced (one species has a stronger effect on the second than the second has on the first
	if (imbalance != 0.5)
		imbalanceCompetition();
	// if competition is transitive (community forms a hierarchy) or intransitive (community has loops breaking hierarchy)
	if (fecundity_transitivity_type != 0 || growth_transitivity_type!=0)
		setCompetitionTransitivity();

	// RNG discard after drawing random parameters, as some use rejection sampling
	unsigned long long max_random_count = (unsigned long long) 1000. * (4. * num_species * num_species );
	if (random_count > max_random_count) {
		if (id == 0)
			fprintf(stderr, "Error, too many random numbers used to generate competition and parameters.\n");
			fprintf(stderr, "Probable causes are the parameterization of TNormal distribution or severe competition correlation\n");
		MPI_Finalize();
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

	// read in delta (competition distance) from parameter file and all other parameters from competition file
	getParameter(delta, num_species, "Delta", 2);
	loadCompetition();

	// RNG discard after loading random parameters, as some used rejection sampling
	unsigned long long max_random_count = (unsigned long long) 1000. * (4. * num_species * num_species );
	if (random_count > max_random_count) {
		if (id == 0)
			fprintf(stderr, "Error, too many random numbers used to generate competition and parameters.\n");
			fprintf(stderr, "Probable causes are the parameterization of TNormal distribution or severe competition correlation\n");
		MPI_Finalize();
		exit(0);
	}
	discardRandom(max_random_count - random_count);

	initializeLattice();

	return;
}

void Simulation::initializeRestartSimulation() {
	/* method used if a previous simulation failed before completing or if you want to extend the simulation.
	continues the simulation from where it left off. this method initializes the lattice, reloads the parameters from
	the previous simulation, and starts up the RNG for the appropriate time step. restarted simulations are identical
	in both parameters and dynamics. */

	// load seeds from competition file and sends to RNG
	loadSeeds();
	seedGenerator();
	// load dispersal from file
	loadDispersal();
	// read in delta (competition distance) from parameter file and all other parameters from competition file
	getParameter(delta, num_species, "Delta", 2);
	loadCompetition();
	loadLattice();

	return;
}

void Simulation::reinitializeSimulation(int time_step) {
	/* method used for two purposes: 1) to reinitialize worker copies of the simulation given the random seeds drawn
	by the central task, and 2) to reinitialize the simulation after the number of species drops below the minimum persistence.
	method deletes already written output files (only for purpose 2) */

	int i;

	// delete previously written output files if simulation is being reinitialized after persistence is too low
	if (id == 0) {
		std::string fname;
			fname = outfile_base + "_competition.csv";
			remove(fname.c_str());
		for (i = 0; i < time_step; i++) {
			fname = outfile_base+"_" + std::to_string(i) + ".csv";
			remove(fname.c_str());
		}
		for (i = 1; i < num_species + 1; i++) { 
			fname = outfile_base + "_dispersal_s" + std::to_string(i) + "_0.csv";
			remove(fname.c_str());
			fname = outfile_base + "_dispersal_s" + std::to_string(i) + "_1.csv";
			remove(fname.c_str());
		}
		num_restarts++;
	}
	random_count = 0;

	if (competition_filename.size() == 0 || restart_time != 0)
		// new simulation, seeds from id = 0
		initializeRandomSimulation();
	else if (restart_time == 0)
		// replicate simulation, seeds from id = 0
		initializeReplicateSimulation();

	// calculate properties of pairwise and community-level interactions
	getSpeciesAbundance();
	getImbalanceMean();
	getDiscreteTransitivity();
	getFecundityGrowthCorrelation();

	// RNG is set to jump ahead to the maximum number of random draws that could have already been used by the RNG when initializing the simulation
	unsigned long long max_random_count = (unsigned long long) 1000. * (4. * num_species * num_species + 5. * lattice_size * lattice_size);
	if (random_count > max_random_count) {
		if (id == 0)
			fprintf(stderr, "Error, too many random numbers used to generate initial conditions.\n");
			fprintf(stderr, "Probable causes are the parameterization of TNormal distribution or severe competition correlation\n");
		MPI_Finalize();
		exit(0);
	}

	discardRandom(max_random_count - random_count);
	random_count = 0;

	return;
}

void Simulation::allocSimulation() {
	/* allocate memory for all arrays used in simulations, including parameter arrays, lattice with species locations, and dispersal lattice with seed locations. */

	int i, j, k;
	species_abundance = new int[num_species];
	delta = new int[num_species];
	juvenile_survival_probability = new double[num_species];
	adult_survival_probability = new double[num_species];
	maximum_competition = new double[num_species];
	dispersal_probability = new double[num_species];
	if (!species_abundance || !delta || !juvenile_survival_probability || !adult_survival_probability || !maximum_competition) {
		fprintf(stderr, "Error, unable to allocate memory for survival or maximum competition arrays\n");
		MPI_Finalize();
		exit(-1);
	}
	fecundity_row_sum = new double[num_species];
	growth_row_sum = new double[num_species];
	fecundity_transitivity = new double *[num_species];
	growth_transitivity = new double *[num_species];
	if (!fecundity_row_sum || !growth_row_sum || !fecundity_transitivity || !growth_transitivity) {
		fprintf(stderr, "Error, unable to allocate memory for transitivity arrays\n");
		MPI_Finalize();
		exit(-1);
	}
	species_occupancy = new double[num_species];
	dispersal_length = new double[num_species];
	intrinsic_fecundity = new double[num_species];
	if (!species_occupancy  || !dispersal_length|| !intrinsic_fecundity ) {
		fprintf(stderr, "Error, unable to allocate memory for species occupancy, dispersal, or intrinsic fecundity arrays\n");
		MPI_Finalize();
		exit(-1);
	}
	competition_fecundity = new double *[num_species];
	competition_growth = new double *[num_species];
	if (!competition_fecundity || !competition_growth) {
		fprintf(stderr, "Error, unable to allocate memory for growth or fecundity competition matrices\n");
		MPI_Finalize();
		exit(-1);
	}
	for (i = 0; i < num_species ; i++) {
		species_abundance[i] = 0;
		delta[i] = 0;
		juvenile_survival_probability[i] = 0.;
		adult_survival_probability[i] = 0.;
		maximum_competition[i] = 0.;
		dispersal_probability[i] = 0.;
		species_occupancy[i] = 0.;
		dispersal_length[i] = 0.;
		intrinsic_fecundity[i] = 0.;
		fecundity_row_sum[i] = 0.;
		growth_row_sum[i] = 0.;

		competition_fecundity[i] = new double[num_species];
		competition_growth[i] = new double[num_species];
		fecundity_transitivity[i] = new double[num_species];
		growth_transitivity[i] = new double[num_species];
		if (!fecundity_transitivity[i] || !growth_transitivity[i]) {
			fprintf(stderr, "Error, unable to allocate memory for transitivity arrays\n");
			MPI_Finalize();
			exit(-1);
		}
		if (!competition_fecundity[i] || !competition_growth[i]) {
			fprintf(stderr, "Error, unable to allocate memory for competition matrices\n");
			MPI_Finalize();
			exit(-1);
		}
		for (j = 0; j < num_species; j++) {
			fecundity_transitivity[i][j] = 0.;
			growth_transitivity[i][j] = 0.;
			competition_fecundity[i][j] = 0.;
			competition_growth[i][j] = 0.;
		}
	}

	lattice = new int *[lattice_size];
	next_lattice = new int *[lattice_size];
	dispersal_lattice = new double **[lattice_size];
	next_dispersal_lattice = new double **[lattice_size];
	if (!lattice || !next_lattice || !dispersal_lattice || !next_dispersal_lattice) {
		fprintf(stderr, "Error, unable to allocate memory for lattice and dispersal lattice\n");
		MPI_Finalize();
		exit(-1);
	}
	for (i = 0; i < lattice_size; i++) {
		lattice[i] = new int[lattice_size];
		next_lattice[i] = new int[lattice_size];
		dispersal_lattice[i] = new double *[lattice_size];
		next_dispersal_lattice[i] = new double *[lattice_size];
		if (!lattice[i] || !next_lattice[i] || !dispersal_lattice[i] || !next_dispersal_lattice[i]) {
			fprintf(stderr, "Error, unable to allocate memory for lattice and dispersal lattice\n");
			MPI_Finalize();
			exit(-1);
		}
		for (j = 0; j < lattice_size; j++) {
			lattice[i][j] = 0;
			next_lattice[i][j] = 0;
			dispersal_lattice[i][j] = new double[num_species];
			next_dispersal_lattice[i][j] = new double[num_species];
			if (!dispersal_lattice[i][j] || !next_dispersal_lattice[i][j]) {
				fprintf(stderr, "Error, unable to allocate memory for dispersal lattice\n");
				MPI_Finalize();
				exit(-1);
			}
			for (k = 0; k < num_species; k++) {
				dispersal_lattice[i][j][k] = 0;
				next_dispersal_lattice[i][j][k] = 0;
			}
		}
	}
	return;
}


void Simulation::initializeLattice() {
	/* initializes lattice with randomly located species, depending on the occupancy probability, which defines both the total
	occupancy of the lattice and the species specific probabilities. also determines whether individuals are juveniles or adults
	with equal probability. */

	int i, j;

	std::bernoulli_distribution stage_dist(0.5);
	std::bernoulli_distribution occupy_dist(initial_occupancy);
	std::discrete_distribution<int> species_dist(species_occupancy, species_occupancy + num_species);

	for (i = 0; i < lattice_size; i++) {
		for (j = 0; j < lattice_size; j++) {
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
	/* calculate the number of individuals of each species in the lattice. */

	int i, j, s;

	for (i = 0; i < num_species; i++)
		species_abundance[i] = 0;

	for (i = 0; i < lattice_size; i++) {
		for (j = 0; j < lattice_size; j++) {
			s = abs(lattice[i][j]);
			if (s != 0)
				species_abundance[s - 1]++;
		}
	}
	return;
}


void Simulation::resetLattice() {
	/* for the current time step, set all elements of the lattice and the dispersal lattice to 0. used in simulation for the first time step (t = 0). */

	int i, j, k;

	for (i = 0; i < lattice_size; i++) {
		for (j = 0; j < lattice_size; j++) {
			lattice[i][j] = 0;
			for (k = 0; k < num_species; k++) {
				dispersal_lattice[i][j][k] = 0.;
			}
		}
	}
	return;
}


void Simulation::resetNextLattice() {
	/* for the next time step, set all elements of the lattice and the dispersal lattice to 0. used in simulation so that workers
	can reset their local copies of 'next_lattice' and 'next_dispersal_lattice' */

	int i, j, k;

	for (i = 0; i < lattice_size; i++) {
		for (j = 0; j < lattice_size; j++) {
			next_lattice[i][j] = 0;
			for (k = 0; k < num_species; k++) {
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
}



void Simulation::updateSingleSite(int i, int j) {
	/* this method runs through all processes, including germination, survival, growth, reproduction, and death,
	for a single site in the lattice. workers use this method to update their local copies of 'next_lattice' and 'next_dispersal_lattice' */
	std::chrono::steady_clock::time_point t_start = std::chrono::steady_clock::now();

	int k, l, kp, lp, pathway = 0;
	unsigned long long start_random_count = random_count;

	int this_delta = 0;
	int this_species = 0;
	int this_stage = 0;
	int this_max_dispersal = 0;
	double this_survival_probability = 0;
	double this_intrinsic_fecundity = 0.;
	double this_dispersal_probability = 0.;
	double this_dispersal_length = 0.;
	double this_maximum_competition = 0.;

	double total_seeds = 0;
	int total_abundance = 0;

	double growth_probability = 0;
	double this_fecundity = 0.;
	double distance_probability_sum = 0.;

	int *neighborhood;
	int *neighborhood_abundance;
	double **distance_probability;

	this_species = abs(lattice[i][j]);  // species in this site i, j
	
	// if the site is open, determine whether or not something germinates with a Bernoulli probability based on the total number of seeds in the site
	// if something germinates, select the species that germinates with probability equal to the relative abundance of each species's seeds in this site.
	if (this_species == 0) {
		total_seeds = 0;
	
		for (k = 0; k < num_species; k++)
			total_seeds += dispersal_lattice[i][j][k];

		std::bernoulli_distribution germ_dist(1. - pow(1.- germination_probability, total_seeds));

		if (germ_dist(generateRandom())) {
			pathway++;
			std::discrete_distribution<int> species_dist(dispersal_lattice[i][j], dispersal_lattice[i][j] + num_species);
			l = abs(species_dist(generateRandom()));
			next_lattice[i][j] = -(l + 1);
			species_abundance[l]++;
		}
		else {
			next_lattice[i][j] = 0;
		}
	}
	else {
		pathway++;
		// if the site is not open, determine whether it's occupied by a juvenile or an adult
		this_stage = lattice[i][j] / abs(lattice[i][j]);
		// the individual survives with stage-specific probability
		if (this_stage > 0)
			this_survival_probability = adult_survival_probability[this_species - 1] ;
		else
			this_survival_probability = juvenile_survival_probability[this_species - 1];
		// determine intrinsic fecundity, dispersal probability, dispersal length, and maximum competition based on the species identity
		this_intrinsic_fecundity = intrinsic_fecundity[this_species - 1];
		this_dispersal_probability = dispersal_probability[this_species - 1];
		this_dispersal_length = dispersal_length[this_species - 1];
		this_maximum_competition = maximum_competition[this_species - 1];

		std::bernoulli_distribution survival_dist(this_survival_probability);
		// if species survives, it persists to next time step
		if (survival_dist(generateRandom())) {
			pathway++;
			next_lattice[i][j] = lattice[i][j];

			// determine abundance of each species in the neighborhood of the focal individual
			// neighborhood limits depend on the parameter delta for this_species
			this_delta = delta[this_species - 1];
			neighborhood = new int[(2 * this_delta + 1) * (2 * this_delta + 1)]; 
			if (!neighborhood) {
					fprintf(stderr, "Error, neighborhood memory allocation failed for site %d %d\n", i, j);	
					MPI_Finalize();
					exit(-1);
			}
			for (k = 0; k < (2 * this_delta + 1) * (2 * this_delta + 1); k++)
				neighborhood[k] = 0;
			neighborhood_abundance = new int[num_species];
			if (!neighborhood_abundance) {
					fprintf(stderr, "Error, neighborhood abundance memory allocation failed for site %d %d\n", i, j);	
					MPI_Finalize();
					exit(-1);
			}
			for (k = 0; k < num_species; k++)
				neighborhood_abundance[k] = 0;	
			for (k = 0; k < 2 * this_delta + 1; k++) {
				for (l = 0; l < 2 * this_delta + 1; l++) {
					if (k != this_delta || l != this_delta) {
						kp = k;
						lp = l;
						if (i - this_delta + kp < 0)
							kp += lattice_size - i;
						else if (i - this_delta + kp > lattice_size - 1)
							kp += -lattice_size;	
						if (j - this_delta + lp < 0)
							lp += lattice_size - j;
						else if (j - this_delta + lp > lattice_size - 1)
							lp += -lattice_size;
						neighborhood[l + (2 * this_delta + 1) * k] = lattice[i - this_delta + kp][j - this_delta + lp];
						if (neighborhood[l + (2 * this_delta + 1) * k] != 0) {
							neighborhood_abundance[abs(neighborhood[l + (2 * this_delta + 1) * k]) - 1]++;
							total_abundance++;
						}	
					}
				}
			}
			delete[] neighborhood;

			// if focal individual is a juvenile, it will grow to become an adult with a probability dependent on competition with individuals in its neighborhood
			// depends on the growth competition matrix, which dictates these interactions. the maximum probability is set as a parameter (never a 100% chance of growing)
			if (this_stage < 0) {
				pathway++;
				for (k = 0; k < num_species; k++)
					growth_probability += competition_growth[this_species - 1][k] * ((double) neighborhood_abundance[k]) / ((double) total_abundance);
				growth_probability = this_maximum_competition * exp(growth_probability);
				if (growth_probability > this_maximum_competition)
					growth_probability = this_maximum_competition;

				std::bernoulli_distribution stage_dist(growth_probability);

				if (stage_dist(generateRandom()))
					next_lattice[i][j] = abs(next_lattice[i][j]);
			}
			else {
				pathway+=2;
				// if focal individual is an adult, it will reproduce with a fecundity based on competition with individuals in its neighborhood
				if (total_abundance != 0) {
					for (k = 0; k < num_species; k++)
						this_fecundity += competition_fecundity[this_species - 1][k] * ((double) neighborhood_abundance[k]) / ((double) total_abundance);
					this_fecundity = this_intrinsic_fecundity * exp(this_fecundity);
				}
				else {
					this_fecundity = this_intrinsic_fecundity;
				}

				// fecundity, the number of seeds produced by the focal individuals, is spread over the entire matrix
				// the number of seeds at each site is an exponential function of Euclidean distance
				// probability of dispersal quickly decays with distance depending on the dispersal length parameter
				// dispersal can occur within two times the dispersal length of the species
				// setting a dispersal limit reduces size of distance_probability array and improves speed
				this_max_dispersal = (int) round(2 * dispersal_length[this_species - 1]);
				distance_probability = new double*[this_max_dispersal + 1];
				if (!distance_probability) {
					fprintf(stderr, "Error, distance probability memory allocation failed for site %d %d\n", i, j);	
					MPI_Finalize();
					exit(-1);
				}
				for (k = 0; k <= this_max_dispersal; k++) {
					distance_probability[k] = new double[this_max_dispersal];
					if (!distance_probability[k]) {
						fprintf(stderr, "Error, distance probability memory allocation failed for site %d %d\n", i, j);	
						MPI_Finalize();
						exit(-1);
					}
					for (l = 1; l <= this_max_dispersal; l++) {
						double distance_x = (double) k;
						double distance_y = (double) l;
						double distance = sqrt(distance_x * distance_x + distance_y * distance_y);
						double probability = 0.;
						if (round(distance) <= this_max_dispersal) {
							probability = exp(log(this_dispersal_probability) / dispersal_length[this_species - 1] * distance);
							distance_probability_sum +=	4. * probability;
						}				
						distance_probability[k][l - 1] = probability;
					}
				}
				for (k = 0; k <= this_max_dispersal; k++) {
					for (l = 1; l <= this_max_dispersal; l++) {
						if (k == 0) {
							int i1 = (i + l) % lattice_size;
							int i2 = (i - l) % lattice_size;
							if (i2 < 0)
								i2 += lattice_size;
							int j1 = (j + l) % lattice_size;
							int j2 = (j - l) % lattice_size;
							if (j2 < 0)
								j2 += lattice_size;
							next_dispersal_lattice[i][j1][this_species - 1] += this_fecundity * distance_probability[k][l - 1] / distance_probability_sum;
							next_dispersal_lattice[i][j2][this_species - 1] += this_fecundity * distance_probability[k][l - 1] / distance_probability_sum;
							next_dispersal_lattice[i1][j][this_species - 1] += this_fecundity * distance_probability[k][l - 1] / distance_probability_sum;
							next_dispersal_lattice[i2][j][this_species - 1] += this_fecundity * distance_probability[k][l - 1] / distance_probability_sum;
							}
						else {
							int i1 = (i + k) % lattice_size;
							int i2 = (i - k) % lattice_size;
							if(i2 < 0)
								i2 += lattice_size;
							int j1 = (j + l) % lattice_size;
							int j2 = (j - l) % lattice_size;
							if(j2 < 0)
								j2 += lattice_size;
							next_dispersal_lattice[i1][j1][this_species - 1] += this_fecundity * distance_probability[k][l - 1] / distance_probability_sum;
							next_dispersal_lattice[i2][j1][this_species - 1] += this_fecundity * distance_probability[k][l - 1] / distance_probability_sum;
							next_dispersal_lattice[i1][j2][this_species - 1] += this_fecundity * distance_probability[k][l - 1] / distance_probability_sum;
							next_dispersal_lattice[i2][j2][this_species - 1] += this_fecundity * distance_probability[k][l - 1] / distance_probability_sum;
						}
					}
					delete[] distance_probability[k];
				}
				delete[] distance_probability;
			}
			delete[] neighborhood_abundance;
		}
		else {
			pathway++;
			// focal individual dies
			next_lattice[i][j] = 0;
			species_abundance[this_species - 1]--;
		}
	}
	// no matter what happened in this site, a total of four random numbers will be discarded (the maximum number of random numbers used in the simulation)
	discardRandom(((unsigned long long) 4 - (random_count - start_random_count)));

	std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();
	std::chrono::duration<double, std::milli> duration = (t_end - t_start) / 1000.;
	/* if (i < 10 && j < 10) {
		fprintf(stdout, "    Done site %d, %d. Pathway %d in %.4e seconds\n", i, j, pathway, duration);
	} */
	return;
}


void Simulation::nextToThis() {
	/* this method updates the species in the lattice at this time step to the next time step,
	and sets the next time step to be unoccupied. does the same for seeds in the dispersal
	lattice */

	int i, j, k;

	for (i = 0; i < lattice_size; i++) {
		for (j = 0; j < lattice_size; j++) {
			lattice[i][j] = next_lattice[i][j];
			next_lattice[i][j] = 0;
			for (k = 0; k < num_species; k++) {
				dispersal_lattice[i][j][k] = next_dispersal_lattice[i][j][k];
				next_dispersal_lattice[i][j][k] = 0;
			}
		}
	}
	return;
}


void Simulation::saveLattice(int time_step) {
	/* saves the species locations in the lattice from the current time step to file */

	int i, j;	

	std::ofstream lattice_file;
	lattice_file.open(outfile_base+"_" + std::to_string(time_step) + ".csv", std::ios::out | std::ios::trunc);

	if (!lattice_file.is_open()) {
			if (id == 0)
				fprintf(stderr, "Error, could not open time step %d lattice file file to load\n", time_step);
			MPI_Finalize();
			exit(0);
	}
	for (i = 0; i < lattice_size; i++) {
		for (j = 0; j < lattice_size; j++) {
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
			MPI_Finalize();
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

	competition_file << "# Fecundity transitivity row sum:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << fecundity_row_sum[i];
		if (i < num_species - 1)
			competition_file << ",";
	}
	competition_file << std::endl;

	competition_file << "# Growth row sum:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << growth_row_sum[i];
		if (i < num_species - 1)
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


int Simulation::getSpecies() {
	return num_species;
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

	int i, p = 0;

	for (i = 0; i < num_species; i++) {
		if (species_abundance[i] != 0)
			p++;
	}
	return p;
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
