#ifndef SIMULATION_H
#define SIMULATION_H

#include <string>
#include <random>
#include <fstream>
#include <sstream>
#include <cmath>
#include <sys/stat.h>
#include <chrono>
#include <omp.h>
#include <vector>

class Simulation {

	public:
		Simulation(std::string filename, int p_id);
		void run(void);
		int getLatticeSize(void);
		int getNumSpecies(void);
		void printSpeciesAbundance(void);
		void printNextLattice(void);
		double getGerminationProbability(void);
		double getAdultSurvivalProbability(int s);
		double getJuvenileSurvivalProbability(int s);
		double getCompetitionGrowth(int s1, int s2);
		double getCompetitionFecundity(int s1, int s2);
		double getDelta(int s);
		double getDispersalLength(int s);
		double getDispersalProbability(int s);
		double getIntrinsicFecundity(int s);
		double getMaximumCompetition(int s);
		int getMaxTimeStep(void);
		int getSite(int i, int j);
		int getNextSite(int i, int j);
		double getDispersal(int i, int j, int s);
		double getNextDispersal(int i, int j, int s);
		int getContinueTime(void);
		unsigned int getSeed(int i);

		void setSite(int i, int j, int s);
		void setNextSite(int i, int j, int s);
		void incrementSpeciesAbundance(int s);
		void decrementSpeciesAbundance(int s);
		void addSite(int i, int j, int s);
		void setDispersal(int i, int j, int s, double t);
		void addDispersal(int i, int j, int s, double t);
		void addNextDispersal(int i, int j, int s, double t);
		void resetLattice(void);
		void resetNextLattice(void);
		void setSeed(int t_i, unsigned int t_s);

		void discardRandom(unsigned long long t_num_discard);
		unsigned long long getRandomCount(void);
		void saveLattice(int t_time_step);
		void saveCompetition(void);
		void saveDispersal(int t_time_step);
		void nextToThis(void);
		double getRandomUniformReal(double t_lower_bound, double t_upper_bound);
		double getRandomNormal(double t_mean, double t_sdev);
		int getPersistence(void);
		int getMinPersistence(void);
		void setRandomSeeds(void);
		void reinitializeSimulation(int t_time_step);

		void seedGenerator(std::mt19937& t_random_generator);
		std::mt19937& generateRandom(unsigned long long& t_random_count, std::mt19937& t_random_generator);
		void discardRandom(unsigned long long t_num_discard, unsigned long long & t_random_count, std::mt19937& t_random_generator);
		unsigned long long getMaxRandomCount(void);

	private:

		// Lattice, dispersal, and time step parameters
		unsigned int seeds[5];
		std::vector<int> species_abundance;
		std::mt19937 global_random_generator;
		std::string parameter_filename, outfile_base, outfile_dir;
		int id, num_species, lattice_size, num_steps, max_time_step, continue_time, min_persistence, num_restarts, num_threads;
		unsigned long long random_count, max_random_count;
		double germination_probability, initial_occupancy;
		std::vector<std::vector<int>> lattice, next_lattice;
		std::vector<std::vector<std::vector<double>>> dispersal_lattice, next_dispersal_lattice;

		// Species specific parameters
		std::vector<int> delta, neighborhood_size;
		std::vector<double> species_occupancy, juvenile_survival_probability, adult_survival_probability;
		std::vector<double> maximum_competition, dispersal_probability, dispersal_length, intrinsic_fecundity;
	
		// Competition parameters
		std::string competition_filename, competition_type;
		double competition_lower_bound, competition_diag_lower_bound, competition_upper_bound, competition_diag_upper_bound; 
		double competition_mean, competition_diag_mean, competition_sdev, competition_diag_sdev;
		double competition_correlation, imbalance, fecundity_transitivity_type, growth_transitivity_type, fecundity_growth_relative_hierarchy;
		double fecundity_imbalance_mean, growth_imbalance_mean, fecundity_growth_correlation;
		double fecundity_relative_intransitivity, growth_relative_intransitivity;
		std::vector<std::vector<double>> competition_fecundity, competition_growth, fecundity_transitivity, growth_transitivity; 
		std::vector<double> fecundity_rank, growth_rank;

		void initializeRandomSimulation(void);
		void initializeReplicateSimulation(void);
		void initializeContinueSimulation(void);
		void loadSeeds(void);
		void loadLattice(void);
		void loadDispersal(void);
		void loadCompetition(void);

		// Allocation and random seed
		std::mt19937& generateRandom(void);
		unsigned int getRandom(void);
		void allocSimulation(void);
		void initializeLattice(void);

		// Read input and set random parameters
		void checkInputFormat(void);
		void seedGenerator(void);
		int getParameter(int & t_value, std::string t_parameter_name, int t_essential);
		void getParameter(double & t_value, std::string t_parameter_name, int t_essential);
		void getParameter(std::string & t_value, std::string t_parameter_name, int t_essential);
		void getParameter(std::vector<int> & t_vector, std::string t_parameter_name, int t_essential);
		void getParameter(std::vector<double> & t_vector, std::string t_parameter_name, int t_essential);
		void initializeNormalRandomArray(std::vector<double> & t_vector, std::vector<double> & t_mean, std::vector<double> & t_sdev);
		void setRandomParameter(std::vector<double> & t_vector, std::string t_parameter_name, int t_type);
		std::string trimString(std::string t_string);
		std::string trimStringNoComment(std::string t_string);

		// Competition Initialization 
		void initializeUniformCompetition(void);
		void initializeTNormalCompetition(void);
		void initializeUniformCorrelatedCompetition(void);
		void initializeTNormalCorrelatedCompetition(void);
		void imbalanceCompetition(void);
		void setCompetitionTransitivity(void);
		void setGrowthCompetitionTransitivity(void);
		void shuffleArray(std::vector<int> & t_vector);		
		void shuffleMatrix(std::vector<std::vector<double>> & t_vector);
		// Competition properties
		void checkCorrelation(void);
		void getSpeciesAbundance(void);
		void getFecundityGrowthCorrelation(void);
		void getImbalanceMean(void);
		void getDiscreteTransitivity(void);
		void getDiscreteFecundityTransitivity(void);
		void getDiscreteGrowthTransitivity(void);
		void getContinuousTransitivity(void);

};

#endif
