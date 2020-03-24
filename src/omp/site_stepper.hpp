#include "simulation.hpp"
#include <vector>

class SiteStepper {

	public:
	SiteStepper(Simulation & t_sim, std::mt19937& t_random_generator, unsigned long long& t_random_count);
	void updateSingleSite(int t_i, int t_j);

	private:
	Simulation * sim;
	std::mt19937 * random_generator;
	unsigned long long * random_count;

	std::bernoulli_distribution bernoulli_distribution;
	std::discrete_distribution<int> discrete_distribution;
	std::vector<double> discrete_distribution_weights;
	std::vector<int> deltas;
	std::vector<int> neighborhood_lengths;
	std::vector<int> neighborhood_sizes;
	std::vector<int> neighborhood_abundances;
	std::vector<std::vector<std::vector<double>>> dispersals;

	int i, j, k, l, kp, lp, i1, i2, j1, j2;
	unsigned long long start_random_count;

	int lattice_size;
	int num_species;
	double germination_probability;
	int species = 0;
	int species_index = 0;
	int stage = 0;
	int neighbor = 0;
	double bernoulli_probability = 0.;

	double total_seeds = 0;
	int total_abundance = 0;
	double competition_fecundity_factor = 0.;

	void initializeDispersals(void);
	void updateEmptySite(void);
	bool drawSurvival(void);
	void countNeighborhoodAbundances(void);
	void updateJuvenileSite(void);
	void disperseSeeds(void);

};

