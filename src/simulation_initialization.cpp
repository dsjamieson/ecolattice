
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2020			
	 *
	 *	methods for Ecolattice class. allocate vectors for the 
	 *	simulation and initializes simulation parameters and lattice.
	 *	includes methods for the various types of simulations:
	 *	random (default), continue, repeat, and replicate.
	 *
	 ***********************************************************/

#include "ecolattice.hpp"

void Ecolattice::allocSimulation() {
	/* allocate memory for all vectors used in simulations, including parameter arrays, 
		lattice with species locations, and dispersal lattice with seed locations. */

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

void Ecolattice::initializeRandomSimulation() {
	/* initializes the simulation lattice with species locations, and draws random variates for species-specific parameters
	(dispersal, competition, etc.). also checks that parameter values are appropriate. */

	// send random seeds to RNG
	seedRandomGenerator();

	// set species specific parameters, potentially random
	getParameter(delta, "Delta", 2);
	for (int i = 0; i < num_species; i++) {
		neighborhood_size[i] = (2 * delta[i] + 1) * (2 * delta[i] + 1) - 1;
	}
	initializeRandomParameter(species_occupancy, "SpeciesOccupancy", 3);
	initializeRandomParameter(juvenile_survival_probability, "JuvenileSurvival", 2);
	initializeRandomParameter(adult_survival_probability, "AdultSurvival", 2);
	initializeRandomParameter(maximum_competition, "MaximumCompetition", 2);
	initializeRandomParameter(dispersal_probability, "DispersalProbability", 2);
	initializeRandomParameter(dispersal_length, "DispersalLength", 4);
	for (int i = 0; i < num_species; i++) {
		if ((int) round(2 * dispersal_length[i]) > lattice_size) {
			if (id == 0)
				fprintf(stderr, "Error, DispersalLength must be less than half of LatticeSize\n");
			exit(0);
		}
	}
	initializeRandomParameter(intrinsic_fecundity, "Fecundity", 4);
	
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

void Ecolattice::initializeReplicateSimulation() {
	/* method used to start a new simulation with the same random parameters, used for replicates.
	uses competition matrices, fecundities, occupancies, etc. from file, specified in 
	"CompetitionFile." replicates have the same parameters, but are different realizations
	(i.e., dynamics will differ). */

	// send seeds to RNG
	seedRandomGenerator();
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

void Ecolattice::initializeContinueSimulation() {
	/* method used if a previous simulation failed before completing or if you want to extend the simulation.
	continues the simulation from where it left off. this method initializes the lattice, reloads the parameters from
	the previous simulation, and starts up the RNG for the appropriate time step. restarted simulations are identical
	in both parameters and dynamics. */

	// load seeds from competition file and sends to RNG
	loadSeeds();
	seedRandomGenerator();
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

void Ecolattice::initializeLattice() {
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
