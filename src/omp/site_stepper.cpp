
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2018			
	 *
	 *		methods for the Simulation class, for omp version
	 *		of the code.
	 *
	 ***********************************************************/

#include "site_stepper.hpp"

SiteStepper::SiteStepper(Simulation & t_sim, std::mt19937 & t_random_generator, unsigned long long & t_random_count) {
	sim = &t_sim;
	random_generator = &t_random_generator;
	random_count = &t_random_count;
	num_species = sim->getNumSpecies();
	lattice_size = sim->getLatticeSize();
	germination_probability = sim->getGerminationProbability();
	discrete_distribution_weights.resize(num_species);
	deltas.resize(num_species);
	neighborhood_lengths.resize(num_species);
	neighborhood_sizes.resize(num_species);
	neighborhood_abundances.resize(num_species);
	for(k = 0; k < num_species; k++) {
		deltas[k] = sim->getDelta(k);
		neighborhood_lengths[k] = 2 * deltas[k] + 1;
		neighborhood_sizes[k] = neighborhood_lengths[k] * neighborhood_lengths[k] - 1;
	} 
	initializeDispersals();
}

void SiteStepper::initializeDispersals() {

	// fecundity, the number of seeds produced by the focal individuals, is spread over the entire matrix
	// the number of seeds at each site is an exponential function of Euclidean distance
	// probability of dispersal quickly decays with distance depending on the dispersal length parameter
	// dispersal can occur within two times the dispersal length of the species
	// setting a dispersal limit reduces size of distance_probability array and improves speed		



	double distance;
	std::vector<double> dispersal_sums(num_species, 0.);
	dispersals.resize(num_species);
	for(k = 0; k < num_species; k++) {
		dispersals[k].resize((int) round(2*sim->getDispersalLength(k) + 1.));
		for(unsigned long i = 0; i < dispersals[k].size(); i++) {
			dispersals[k][i].resize(dispersals[k].size() - 1);
			for(unsigned long j = 0; j < dispersals[k][i].size(); j++) {
				distance = sqrt((double) i * i + (j + 1) * (j + 1));
				if(round(distance) <= dispersals[k].size() - 1) {
					dispersals[k][i][j] = exp(log(sim->getDispersalProbability(k)) * distance / sim->getDispersalLength(k));
					dispersal_sums[k] += 4. * dispersals[k][i][j];
				}
				else {
					dispersals[k][i][j] = 0.;
				}
			}
		}

		for(unsigned long i = 0; i < dispersals[k].size(); i++) {
			for(unsigned long j = 0; j < dispersals[k][i].size(); j++) {
					dispersals[k][i][j] *= sim->getIntrinsicFecundity(k) / dispersal_sums[k];
			}
		}
	}
	return;
}


void SiteStepper::updateSingleSite(int t_i, int t_j) {
	start_random_count = *random_count;
	i = t_i;
	j = t_j;
	//fprintf(stdout, "%d %d: %d\n", i, j, sim->getSite(i, j));
	species = abs(sim->getSite(i, j));
	if (species == 0) {
		updateEmptySite();
	}
	else {
		// if the site is not open, determine whether it's occupied by a juvenile or an adult
		stage = sim->getSite(i, j)/species;
		species_index = species - 1;
		if (drawSurvival()) {
			sim->setNextSite(i, j, sim->getSite(i, j));
			countNeighborhoodAbundances();
			if (stage < 0) {
				updateJuvenileSite();
			}
			else {
				disperseSeeds();
			}
		}
		else {
			sim->decrementSpeciesAbundance(species_index);
		}
	}
	// no matter what happened in this site, a total of four random numbers will be discarded (the maximum number of random numbers used in the simulation)
	sim->discardRandom(((unsigned long long) 4 - (*random_count - start_random_count)), *random_count, *random_generator);
	return;
}

void SiteStepper::updateEmptySite(void) {
	total_seeds = 0.;
	for (k = 0; k < num_species; k++)
		total_seeds += sim->getDispersal(i, j, k);
	bernoulli_probability = 1. - pow(1. - germination_probability, total_seeds);
	bernoulli_distribution.param((std::bernoulli_distribution::param_type) bernoulli_probability);
	if (bernoulli_distribution(sim->generateRandom(*random_count, *random_generator))) {
		// Change when vectors are included for Simulation
		for(k = 0; k < num_species; k++)
			discrete_distribution_weights[k] = sim->getDispersal(i,j,k);
		discrete_distribution.param(std::discrete_distribution<int>::param_type(discrete_distribution_weights.begin(), discrete_distribution_weights.end()));
		l = discrete_distribution(sim->generateRandom(*random_count, *random_generator));
		sim->setNextSite(i, j, -(l + 1));
		sim->incrementSpeciesAbundance(l);
	}
	return;
}

bool SiteStepper::drawSurvival(void) {
	// the individual survives with stage-specific probability
	if (stage > 0)
		bernoulli_probability = sim->getAdultSurvivalProbability(species_index);
	else
		bernoulli_probability = sim->getJuvenileSurvivalProbability(species_index);
	bernoulli_distribution.param((std::bernoulli_distribution::param_type) bernoulli_probability);
	// if species survives, it persists to next time step
	return bernoulli_distribution(sim->generateRandom(*random_count, *random_generator));
}

void SiteStepper::countNeighborhoodAbundances(void) {
	// determine abundance of each species in the neighborhood of the focal individual
	// neighborhood limits depend on the parameter delta for species
	for (k = 0; k < num_species; k++) {
		neighborhood_abundances[k] = 0;
	}
	total_abundance = 0;
	for (k = 0; k < neighborhood_lengths[species_index]; k++) {
		for (l = 0; l < neighborhood_lengths[species_index]; l++) {
			// skip center case
			if (k != deltas[species_index] || l != deltas[species_index]) {
				// lattice indices for site neighbors, accounting for periodic boundary conditions
				kp = k;
				lp = l;
				// if lattice row index is negative, wrap around to end of row
				if (i - deltas[species_index] + kp < 0)
					kp += lattice_size - i;
				// if lattice row index is past the end of the row, wrap around to start of row
				else if (i - deltas[species_index] + kp > lattice_size - 1)
					kp += -lattice_size;
				// if lattice column index is negative, wrap around to end of column
				if (j - deltas[species_index] + lp < 0)
					lp += lattice_size - j;
				// if lattice column index is past the end of the column, wrap around to start of column
				else if (j - deltas[species_index] + lp > lattice_size - 1)
					lp += -lattice_size;
				// use neighbor lattice indices to place neighbors in neighborhood vector
				neighbor = sim->getSite(i - deltas[species_index] + kp, j - deltas[species_index] + lp);
				if (neighbor != 0) {
					neighborhood_abundances[abs(neighbor) - 1]++;
					total_abundance++;
				}	
			}
		}
	}
	return;
}

void SiteStepper::updateJuvenileSite(void) {
	// if focal individual is a juvenile, it will grow to become an adult with a probability dependent on competition with individuals in its neighborhood
	// depends on the growth competition matrix, which dictates these interactions. the maximum probability is set as a parameter (never a 100% chance of growing)
	bernoulli_probability = 0.;
	for (k = 0; k < num_species; k++)
		bernoulli_probability += sim->getCompetitionGrowth(species_index, k) * ((double) neighborhood_abundances[k]) / ((double) neighborhood_sizes[species_index]);
	bernoulli_probability = std::min(exp(bernoulli_probability), sim->getMaximumCompetition(species_index));
	bernoulli_distribution.param((std::bernoulli_distribution::param_type) bernoulli_probability);
	if (bernoulli_distribution(sim->generateRandom(*random_count, *random_generator)))
		sim->setNextSite(i, j, species);
	return;
}

void SiteStepper::disperseSeeds(void) {
	// if focal individual is an adult, it will reproduce with a fecundity based on competition with individuals in its neighborhood
	competition_fecundity_factor = 0.;
	if (total_abundance != 0) {
		for (k = 0; k < num_species; k++)
			competition_fecundity_factor += sim->getCompetitionFecundity(species_index, k) * ((double) neighborhood_abundances[k]) / ((double) neighborhood_sizes[species_index]);
		competition_fecundity_factor = exp(competition_fecundity_factor);
	}
	else {
		competition_fecundity_factor  = 1.;
	}
	// start at lower right rectangular quadrant
	// determine the four lattice indices that correspond with each dispersal distance
	for (k = 0; k < (int) dispersals[species_index].size(); k++) {
		for (l = 0; l < (int) dispersals[species_index][k].size(); l++) {
			if (k == 0) {
				// in same row as site i, j (l = distance from center)
				// the four lattice indices are in the cardinal directions
				int i1 = (i + l + 1) % lattice_size;
				int i2 = i - (l + 1);
				if (i2 < 0)
					i2 += lattice_size;
				int j1 = (j + l + 1) % lattice_size;
				int j2 = j - (l + 1);
				if (j2 < 0)
					j2 += lattice_size;
				// same row, i1 to the right
				sim->addNextDispersal(i, j1, species_index, competition_fecundity_factor*dispersals[species_index][k][l]);
				// same row, i2 to the left
				sim->addNextDispersal(i, j2, species_index, competition_fecundity_factor*dispersals[species_index][k][l]);
				// same column, i1 up
				sim->addNextDispersal(i1, j, species_index, competition_fecundity_factor*dispersals[species_index][k][l]);
				// same column, i2 down
				sim->addNextDispersal(i2, j, species_index, competition_fecundity_factor*dispersals[species_index][k][l]);
				}
			else {
				// square quadrants not in the same row as site i, j
				i1 = (i + k) % lattice_size;
				i2 = i - k;
				if (i2 < 0)
					i2 += lattice_size;
				j1 = (j + l + 1) % lattice_size;
				j2 = j - (l + 1);
				if (j2 < 0)
					j2 += lattice_size;
				// upper right quadrant
				sim->addNextDispersal(i1, j1, species_index, competition_fecundity_factor*dispersals[species_index][k][l]);
				// upper left quadrant
				sim->addNextDispersal(i2, j1, species_index, competition_fecundity_factor*dispersals[species_index][k][l]);
				// lower right quadrant
				sim->addNextDispersal(i1, j2, species_index, competition_fecundity_factor*dispersals[species_index][k][l]);
				// lower left quadrant
				sim->addNextDispersal(i2, j2, species_index, competition_fecundity_factor*dispersals[species_index][k][l]);
			}
		}
	}
	return;
}
