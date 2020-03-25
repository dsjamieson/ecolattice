
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2020			
	 *
	 *		public methods for the Ecolattice class, used in 
	 *		SiteStepper object to access objects in Ecolattice.
	 *
	 ***********************************************************/

#include "ecolattice.hpp"

unsigned long long Ecolattice::getMaxRandomCount(void) {
	return max_random_count;
}

void Ecolattice::seedRandomGenerator(std::mt19937& this_random_generator) {
	/* seeds a thread's RNG with the object's array of seeds. */
    std::seed_seq seq(seeds, seeds + 5);
	std::vector<std::uint32_t> seed_vector(std::mt19937::state_size);
    seq.generate(seed_vector.begin(), seed_vector.end());
	std::seed_seq seq2(seed_vector.begin(), seed_vector.end());
	this_random_generator.seed(seq2);
	return;
}

int Ecolattice::getNumSpecies(void) {
	return num_species;
}

int Ecolattice::getLatticeSize(void) {
	return lattice_size;
}

double Ecolattice::getGerminationProbability(void) {
	return germination_probability;
}

double Ecolattice::getDelta(int s) {
	return delta[s];
}

double Ecolattice::getMaximumCompetition(int s) {
	return maximum_competition[s];
}

double Ecolattice::getDispersalLength(int s) {
	return dispersal_length[s];
}

double Ecolattice::getDispersalProbability(int s) {
	return dispersal_probability[s];
}

double Ecolattice::getIntrinsicFecundity(int s) {
	return intrinsic_fecundity[s];
}

double Ecolattice::getAdultSurvivalProbability(int s) {
	return adult_survival_probability[s];
}

double Ecolattice::getJuvenileSurvivalProbability(int s) {
	return juvenile_survival_probability[s];
}

double Ecolattice::getCompetitionGrowth(int s1, int s2) {
	return competition_growth[s1][s2];
}

double Ecolattice::getCompetitionFecundity(int s1, int s2) {
	return competition_fecundity[s1][s2];
}

int Ecolattice::getSite(int i, int j) {
	return lattice[i][j];
}

double Ecolattice::getDispersal(int i, int j, int s) {
	return dispersal_lattice[i][j][s];
}

void Ecolattice::incrementSpeciesAbundance(int s) {
	#pragma omp atomic
	species_abundance[s]++;
	return;
}

void Ecolattice::decrementSpeciesAbundance(int s) {
	#pragma omp atomic
	species_abundance[s]--;
	return;
}

void Ecolattice::setNextSite(int i, int j, int s) {
	next_lattice[i][j] = s;
	return;
}

void Ecolattice::addNextDispersal(int i, int j, int s, double t) {
	#pragma omp atomic
	next_dispersal_lattice[i][j][s] += t;
	return;
}

