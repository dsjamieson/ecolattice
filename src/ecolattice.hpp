#ifndef ECOLATTICE_H
#define ECOLATTICE_H

#include <string>
#include <random>
#include <fstream>
#include <sstream>
#include <cmath>
#include <sys/stat.h>
#include <chrono>
#include <vector>
#ifdef OMP
	#include <omp.h>
#endif

class Ecolattice {

	public:
		// constructor and running method
		Ecolattice(std::string filename, int p_id);
		void run(void);
		// public methods used by SiteStepper object
		unsigned long long getMaxRandomCount(void);
		void seedRandomGenerator(std::mt19937& t_random_generator);
		int getNumSpecies(void);
		int getLatticeSize(void);
		double getGerminationProbability(void);
		double getDelta(int s);
		double getMaximumCompetition(int s);
		double getDispersalLength(int s);
		double getDispersalProbability(int s);
		double getIntrinsicFecundity(int s);
		double getAdultSurvivalProbability(int s);
		double getJuvenileSurvivalProbability(int s);
		double getCompetitionGrowth(int s1, int s2);
		double getCompetitionFecundity(int s1, int s2);
		int getSite(int i, int j);
		double getDispersal(int i, int j, int s);
		void incrementSpeciesAbundance(int s);
		void decrementSpeciesAbundance(int s);
		void setNextSite(int i, int j, int s);
		void addNextDispersal(int i, int j, int s, double t);

	private:
		// lattice, dispersal, and time step parameters
		unsigned int seeds[5];
		std::vector<int> species_abundance;
		std::mt19937 global_random_generator;
		std::string parameter_filename, outfile_base, outfile_dir;
		int id, num_species, lattice_size, num_steps, max_time_step, continue_time, min_persistence, num_restarts, num_threads, initialization_scheme;
		unsigned long long random_count, max_random_count;
		double germination_probability, initial_occupancy;
		std::vector<std::vector<int>> lattice, next_lattice;
		std::vector<std::vector<std::vector<double>>> dispersal_lattice, next_dispersal_lattice;
		// species specific parameters
		std::vector<int> delta, neighborhood_size;
		std::vector<double> species_occupancy, juvenile_survival_probability, adult_survival_probability;
		std::vector<double> maximum_competition, dispersal_probability, dispersal_length, intrinsic_fecundity;
		// competition parameters
		std::string competition_filename, competition_type;
		double competition_lower_bound, competition_diag_lower_bound, competition_upper_bound, competition_diag_upper_bound; 
		double competition_mean, competition_diag_mean, competition_sdev, competition_diag_sdev;
		double competition_correlation, imbalance, fecundity_transitivity_type, growth_transitivity_type, fecundity_growth_relative_hierarchy;
		double fecundity_imbalance_mean, growth_imbalance_mean, fecundity_growth_correlation;
		double fecundity_relative_intransitivity, growth_relative_intransitivity;
		std::vector<std::vector<double>> competition_fecundity, competition_growth, fecundity_transitivity, growth_transitivity; 
		std::vector<double> fecundity_rank, growth_rank;

		// read input parameters
		void checkInputFormat(void);
		int getParameter(int & t_value, std::string t_parameter_name, int t_essential);
		void getParameter(double & t_value, std::string t_parameter_name, int t_essential);
		void getParameter(std::string & t_value, std::string t_parameter_name, int t_essential);
		void getParameter(std::vector<int> & t_vector, std::string t_parameter_name, int t_essential);
		void getParameter(std::vector<double> & t_vector, std::string t_parameter_name, int t_essential);
		std::string trimString(std::string t_string);
		std::string trimStringNoComment(std::string t_string);
		// random number methods (random_number_generation)
		void drawRandomSeeds(void);
		void seedRandomGenerator(void);
		void discardRandom(unsigned long long t_num_discard);
		void discardRandom(unsigned long long t_num_discard, unsigned long long & t_random_count, std::mt19937 & t_random_generator);
		std::mt19937& generateRandom(void);
		unsigned int getRandom(void);
		double getRandomUniformReal(double t_lower_bound, double t_upper_bound);
		double getRandomNormal(double t_mean, double t_sdev);
		void initializeRandomParameter(std::vector<double> & t_vector, std::string t_parameter_name, int t_type);
		void initializeNormalRandomArray(std::vector<double> & t_vector, std::vector<double> & t_mean, std::vector<double> & t_sdev);
		// allocation and initialization schemes (simulation_initialization)
		void allocSimulation(void);
		void initializeRandomSimulation(void);
		void initializeReplicateSimulation(void);
		void initializeContinueSimulation(void);
		void initializeLattice(void);
		// loading data from files (load_data)
		void loadSeeds(void);
		void loadLattice(void);
		void loadDispersal(void);
		void loadCompetition(void);
		// saving data to files (save_data)
		void saveCompetition(void);
		void saveLattice(int t_time_step);
		void saveDispersal(int t_time_step);
		// competition initialization (competition_initialization)
		void initializeUniformCompetition(void);
		void initializeTNormalCompetition(void);
		void initializeUniformCorrelatedCompetition(void);
		void initializeTNormalCorrelatedCompetition(void);
		void imbalanceCompetition(void);
		void setCompetitionTransitivity(void);
		void setGrowthCompetitionTransitivity(void);
		void shuffleArray(std::vector<int> & t_vector);		
		void shuffleMatrix(std::vector<std::vector<double>> & t_vector);
		// competition properties (competition_analyzation)
		void getSpeciesAbundance(void);
		void printSpeciesAbundance(void);
		void getFecundityGrowthCorrelation(void);
		void getImbalanceMean(void);
		void getDiscreteTransitivity(void);
		void getDiscreteFecundityTransitivity(void);
		void getDiscreteGrowthTransitivity(void);
		void getContinuousTransitivity(void);
		// run simulation (run)
		int getPersistence(void);
		void nextToThis(void);
		void reinitializeSimulation(int t_time_step);
		
		// debug
		void compareSpeciesAbundance(void);

};

#endif
