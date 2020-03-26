
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2020			
	 *
	 *		defines SiteStepper constructor and associated methods.
	 *		SiteStepper goes through simulation model for each
	 *		site on the lattice.
	 *
	 ***********************************************************/

#include "site_stepper.hpp"

SiteStepper::SiteStepper(Ecolattice & t_sim) {
	/* each thread has the pointer to the simulation object, and accesses
 		 this object to make local copies of recurring parameters. */

	// pointer to Ecolattice
	sim = &t_sim;
	loadEcolattice();
}

void SiteStepper::loadEcolattice(void) {
	// RNG is only accesses within SiteStepper object, and is set up here
	max_random_count = sim->getMaxRandomCount();
	initializeRandomGenerator();
	// access simulation object to get important parameters
	num_species = sim->getNumSpecies();
	lattice_size = sim->getLatticeSize();
	germination_probability = sim->getGerminationProbability();
	discrete_distribution_weights.resize(num_species);
	deltas.resize(num_species);
	neighborhood_lengths.resize(num_species);
	neighborhood_sizes.resize(num_species);
	neighborhood_abundances.resize(num_species);
	max_competition.resize(num_species);
	for (k = 0; k < num_species; k++) {
		deltas[k] = sim->getDelta(k);
		neighborhood_lengths[k] = 2 * deltas[k] + 1;
		neighborhood_sizes[k] = neighborhood_lengths[k] * neighborhood_lengths[k] - 1;
		max_competition[k] = sim->getMaximumCompetition(k);
	} 
	initializeDispersals();
	initializeCompetition();

	return;
}

void SiteStepper::initializeRandomGenerator(void) {
	sim->seedRandomGenerator(random_generator);
	random_count = 0;
	discardRandom(max_random_count);
	random_count = 0;
	return;
}

void SiteStepper::initializeDispersals(void) {
	/* dispersal depends on the dispersal limit parameter and the fecundity of the focal species, 
  		and the Euclidean distance. the only part that depends on the current state of the simulation
  		is the fecundity, which depends on the neighborhood of the focal species. thus, everything except
  		the fecundity competition from the neighborhood can be calculated once by each thread and 
  		referenced in the simulation time step. */ 

	double distance;
	std::vector<double> dispersal_sums(num_species, 0.);
	// dimensions of dispersal array = number species * max dispersal distance * max dispersal distance
	dispersals.resize(num_species);
	for (k = 0; k < num_species; k++) {
		// max dispersal distance is twice the dispersal length parameter of the species (rounded, species-specific)
		dispersals[k].resize((int) round(2 * sim->getDispersalLength(k) + 1.));
		for (unsigned long i = 0; i < dispersals[k].size(); i++) {
			dispersals[k][i].resize(dispersals[k].size() - 1);
			for (unsigned long j = 0; j < dispersals[k][i].size(); j++) {
				// Euclidean distance from focal individual
				distance = sqrt((double) i * i + (j + 1) * (j + 1));
				if (round(distance) <= dispersals[k].size() - 1) {
					// dispersal decays quickly away from the focal individual
					dispersals[k][i][j] = exp(log(sim->getDispersalProbability(k)) * distance / sim->getDispersalLength(k));
					// total probability of seeds over the range of maximum dispersal
					dispersal_sums[k] += 4. * dispersals[k][i][j];
				}
				else {
					dispersals[k][i][j] = 0.;
				}
			}
		}
		// species-specific intrinsic fecundity also does not depend on the simulation state
		for (unsigned long i = 0; i < dispersals[k].size(); i++) {
			for (unsigned long j = 0; j < dispersals[k][i].size(); j++) {
					dispersals[k][i][j] *= sim->getIntrinsicFecundity(k) / dispersal_sums[k];
			}
		}
	}
	return;
}

void SiteStepper::initializeCompetition(void) {
	/* access simulation growth and fecundity competition matrices in SiteStepper by 
 		storing a copy within this scope. */

	adult_survival_probability.resize(num_species);
	juvenile_survival_probability.resize(num_species);
	competition_growth.resize(num_species);
	competition_fecundity.resize(num_species);
	for (k = 0; k < num_species; k++) {
		adult_survival_probability[k] = sim->getAdultSurvivalProbability(k);
		juvenile_survival_probability[k] = sim->getJuvenileSurvivalProbability(k);
		competition_growth[k].resize(num_species);
		competition_fecundity[k].resize(num_species);
		for (l = 0; l < num_species; l++) {
			competition_growth[k][l] = sim->getCompetitionGrowth(k, l);
			competition_fecundity[k][l] = sim->getCompetitionFecundity(k, l);
		}
	}
	return;
}

void SiteStepper::updateSingleSite(int t_i, int t_j, int t_time_step) {
	/* the simulation occurs within this method. for a given time step, determine species
 		and stage at site and apply species and stage specific methods. */

	i = t_i;
	j = t_j;
	// start RNG at correct spot given time step and lattice size (so that code can be parallelized, and repeatable)
	discardRandom(4 * lattice_size * lattice_size * ((unsigned long long) t_time_step - 1) + 4 * (j + lattice_size * i) - random_count);
	start_random_count = random_count;
	// access simulation object to determine the species in site i, j
	species = sim->getSite(i, j);
	// if site is empty
	if (species == 0) {
		updateEmptySite();
	}
	else {
		// if the site is not open, determine whether it's occupied by a juvenile or an adult
		stage = abs(species)/species;
		species = abs(species);
		species_index = species - 1;
		// determine whether species survives given survival probabilities
		if (drawSurvival()) {
			// individual survives, store species in lattice for next time step
			sim->setNextSite(i, j, sim->getSite(i, j));
			// adults and juveniles need to know the abundances of all species in their neighborhood
			countNeighborhoodAbundances();
			if (stage < 0) {
				// focal individual is a juvenile
				updateJuvenileSite();
			}
			else {
				// focl individual is an adult
				disperseSeeds();
			}
		}
		else {
			// individual dies
			sim->decrementSpeciesAbundance(species_index);
		}
	}
	// no matter what happened in this site, a total of four random numbers are discarded (the maximum number of random numbers used in a given site)
	discardRandom(((unsigned long long) 4 - (random_count - start_random_count)));
	return;
}

void SiteStepper::updateEmptySite(void) {
	/* method for empty sites. seeds may germinate, and germinated species depends on seeds that have dispersed to site. */

	total_seeds = 0.;
	// determine the total number of seeds of all species at site i, j
	for (k = 0; k < num_species; k++)
		total_seeds += sim->getDispersal(i, j, k);
	// probability of a seed germinating depends on the total number of seeds
	bernoulli_probability = 1. - pow(1. - germination_probability, total_seeds);
	bernoulli_distribution.param((std::bernoulli_distribution::param_type) bernoulli_probability);
	if (bernoulli_distribution(generateRandom())) {
		// seed germinates, species probabilities are the number of seeds present
		for (k = 0; k < num_species; k++)
			discrete_distribution_weights[k] = sim->getDispersal(i,j,k);
		discrete_distribution.param(std::discrete_distribution<int>::param_type(discrete_distribution_weights.begin(), discrete_distribution_weights.end()));
		l = discrete_distribution(generateRandom());
		// juvenile in next time step
		sim->setNextSite(i, j, -(l + 1));
		sim->incrementSpeciesAbundance(l);
	}
	return;
}

bool SiteStepper::drawSurvival(void) {
	/* draw a random bernouilli variate with probability determined by AdultSurvival and JuvenileSurvival. */

	if (stage > 0)
		bernoulli_probability = adult_survival_probability[species_index]; 
	else
		bernoulli_probability = juvenile_survival_probability[species_index];
	bernoulli_distribution.param((std::bernoulli_distribution::param_type) bernoulli_probability);
	return bernoulli_distribution(generateRandom());
}

void SiteStepper::countNeighborhoodAbundances(void) {
	/* method to  determine abundance of each species in the neighborhood of the focal individual. 
 		method used by both adult and juvenile individuals. */

	for (k = 0; k < num_species; k++) {
		neighborhood_abundances[k] = 0;
	}
	total_abundance = 0;
	// neighborhood limits depend on the parameter delta for species
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
	/* method for juvenile individuals. it will grow to become an adult with a probability dependent
 		 on competition with individuals in its neighborhood. competition depends on the growth competition
		 matrix, which dictates these interactions. the maximum probability is the MaximumCompetition parameter (never a 100% chance of growing) */

	bernoulli_probability = 0.;
	for (k = 0; k < num_species; k++)
		bernoulli_probability += competition_growth[species_index][k] * ((double) neighborhood_abundances[k]) / ((double) neighborhood_sizes[species_index]);
	bernoulli_probability = std::min(exp(bernoulli_probability), max_competition[species_index]);
	bernoulli_distribution.param((std::bernoulli_distribution::param_type) bernoulli_probability);
	// juvenile grows to be an adult in lattice at the next time step
	if (bernoulli_distribution(generateRandom()))
		sim->setNextSite(i, j, species);
	return;
}

void SiteStepper::disperseSeeds(void) {
	/* method for adult individuals. adult fecundity (number of seeds dispersed) depends on competition with 
 		individuals in neighborhood. uses pre-computed factors from initializeDispersal to maximize efficiency
		(these computations are not site-specific).  */

	competition_fecundity_factor = 0.;
	if (total_abundance != 0) {
		// compute the part of fecundity that is site specific (depends on neighborhood)
		for (k = 0; k < num_species; k++)
			competition_fecundity_factor += competition_fecundity[species_index][k] * ((double) neighborhood_abundances[k]) / ((double) neighborhood_sizes[species_index]);
		competition_fecundity_factor = exp(competition_fecundity_factor);
	}
	else {
		// no individuals around focal individual, fecundity = intrinsic fecundity
		competition_fecundity_factor  = 1.;
	}
	// use pre-computed dispersal values from initializeDispersal method
	// going from neighborood indices to lattice indices: start at lower right rectangular quadrant
	// determine the four lattice indices that correspond with each dispersal distance
	for (k = 0; k < (int) dispersals[species_index].size(); k++) {
		for (l = 0; l < (int) dispersals[species_index][k].size(); l++) {
			if (k == 0) {
				// in same row as site i, j (l = distance from center)
				// the four lattice indices are in the cardinal directions
				i1 = (i + l + 1) % lattice_size;
				i2 = i - (l + 1);
				if (i2 < 0)
					i2 += lattice_size;
				j1 = (j + l + 1) % lattice_size;
				j2 = j - (l + 1);
				if (j2 < 0)
					j2 += lattice_size;
				// same row, i1 to the right
				sim->addNextDispersal(i, j1, species_index, competition_fecundity_factor * dispersals[species_index][k][l]);
				// same row, i2 to the left
				sim->addNextDispersal(i, j2, species_index, competition_fecundity_factor * dispersals[species_index][k][l]);
				// same column, i1 up
				sim->addNextDispersal(i1, j, species_index, competition_fecundity_factor * dispersals[species_index][k][l]);
				// same column, i2 down
				sim->addNextDispersal(i2, j, species_index, competition_fecundity_factor * dispersals[species_index][k][l]);
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
				sim->addNextDispersal(i1, j1, species_index, competition_fecundity_factor * dispersals[species_index][k][l]);
				// upper left quadrant
				sim->addNextDispersal(i2, j1, species_index, competition_fecundity_factor * dispersals[species_index][k][l]);
				// lower right quadrant
				sim->addNextDispersal(i1, j2, species_index, competition_fecundity_factor * dispersals[species_index][k][l]);
				// lower left quadrant
				sim->addNextDispersal(i2, j2, species_index, competition_fecundity_factor * dispersals[species_index][k][l]);
			}
		}
	}
	return;
}

std::mt19937& SiteStepper::generateRandom(void) {
	/* add to the random count and get a random draw from the RNG. */
	random_count += 2;
	return random_generator;
}

void SiteStepper::discardRandom(unsigned long long t_num_discard) {
	/* add to the random count and discard values from the RNG */
	random_generator.discard(t_num_discard);
	random_count += t_num_discard;
	return;
}
