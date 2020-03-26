#ifndef SITESTEPPER_H
#define SITESTEPPER_H

#include "ecolattice.hpp"

class SiteStepper {

	public:
		// constructor and method for running
		SiteStepper(Ecolattice & t_sim);
		void updateSingleSite(int t_i, int t_j, int t_time_step);
		void loadEcolattice(void);

	private:
		// objects and counter inherited from Ecolattice object using Ecolattice::run
		Ecolattice * sim;
		std::mt19937  random_generator;
		unsigned long long random_count = 0;
		// vectors and distributions
		std::bernoulli_distribution bernoulli_distribution;
		std::discrete_distribution<int> discrete_distribution;
		std::vector<double> discrete_distribution_weights;
		std::vector<int> deltas;
		std::vector<int> neighborhood_lengths;
		std::vector<int> neighborhood_sizes;
		std::vector<int> neighborhood_abundances;
		std::vector<double> max_competition;
		std::vector<double> adult_survival_probability;
		std::vector<double> juvenile_survival_probability;
		std::vector<std::vector<std::vector<double>>> dispersals;
		std::vector<std::vector<double>> competition_growth;
		std::vector<std::vector<double>> competition_fecundity;
		// interators, simulation, and site properties
		int i, j, k, l, kp, lp, i1, i2, j1, j2;
		unsigned long long max_random_count, start_random_count;
		int lattice_size, num_species;
		double germination_probability;
		int species, species_index, stage, neighbor;
		double bernoulli_probability, total_seeds, competition_fecundity_factor;
		int total_abundance;

		// initialization methods
		void initializeDispersals(void);
		void initializeCompetition(void);
		// time-step methods
		void updateEmptySite(void);
		bool drawSurvival(void);
		void countNeighborhoodAbundances(void);
		void updateJuvenileSite(void);
		void disperseSeeds(void);
		// random number generation methods
		void initializeRandomGenerator(void);
		std::mt19937& generateRandom(void);
		void discardRandom(unsigned long long t_num_discard);

};

#endif
