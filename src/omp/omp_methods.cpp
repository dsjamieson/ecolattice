
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2018			
	 *
	 *		methods for the Simulation class, for omp version
	 *		of the code.
	 *
	 ***********************************************************/

#include "simulation.hpp"


void Simulation::seedGenerator(std::mt19937& this_random_generator) {
	/* create a random vector of seeds (a seed sequence) given the seeds specified randomly or in
	the parameter file (for restarted simulations). seeds fed to the global RNG. */
    std::seed_seq seq(seeds, seeds + 5);
	std::vector<std::uint32_t> seed_vector(std::mt19937::state_size);
    seq.generate(seed_vector.begin(), seed_vector.end());
	std::seed_seq seq2(seed_vector.begin(), seed_vector.end());
	this_random_generator.seed(seq2);
	return;
}

std::mt19937& Simulation::generateRandom(unsigned long long& this_random_count, std::mt19937& this_random_generator) {
	/* add to the random count and get a random draw from the global RNG. */

	this_random_count += 2;
	return this_random_generator;
}

void Simulation::discardRandom(unsigned long long n, unsigned long long& this_random_count, std::mt19937& this_random_generator) {
	/* add to the random count and discard values from the global RNG. */

	this_random_generator.discard(n);
	this_random_count += n;
}

unsigned long long Simulation::getMaxRandomCount() {
	return max_random_count;
}

void Simulation::incrementSpeciesAbundance(int s) {
	#pragma omp atomic
	species_abundance[s]++;
	return;
}

void Simulation::decrementSpeciesAbundance(int s) {
	#pragma omp atomic
	species_abundance[s]--;
	return;
}

void Simulation::addNextDispersal(int i, int j, int s, double t) {
	#pragma omp atomic
	next_dispersal_lattice[i][j][s] += t;
	return;
}

