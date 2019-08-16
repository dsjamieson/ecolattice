#ifndef SIMULATION_H
#define SIMULATION_H

#include <string>
#include <random>
#include <fstream>
#include <sstream>
#include <cmath>
#include <sys/stat.h>
#include <chrono>
#include "simulation.h"
#include <omp.h>

class Simulation {

	public:
		Simulation(std::string filename, int p_id);
		int getLatticeSize();
		int getSpecies();
		int getMaxTimeStep();
		int getSite(int i, int j);
		int getNextSite(int i, int j);
		double getDispersal(int i, int j, int s);
		double getNextDispersal(int i, int j, int s);
		int getRestartTime();
		unsigned int getSeed(int i);

		void setSite(int i, int j, int s);
		void addSite(int i, int j, int s);
		void setDispersal(int i, int j, int s, double t);
		void addDispersal(int i, int j, int s, double t);
		void resetLattice();
		void resetNextLattice();
		void setSeed(int i, unsigned int s);

		void updateSingleSite(int i, int j);
		void discardRandom(unsigned long long n);
		unsigned long long getRandomCount();
		void saveLattice(int time_step);
		void saveCompetition();
		void saveDispersal(int time_step);
		void nextToThis();
		double getRandomUniformReal(double lower_bound, double upper_bound);
		double getRandomNormal(double mean, double sdev);
		int getPersistence();
		int getMinPersistence();
		void setRandomSeeds();
		void reinitializeSimulation(int time_step);

		void seedGenerator(std::mt19937& this_random_generator);
		std::mt19937& generateRandom(unsigned long long& this_random_count, std::mt19937& this_random_generator);
		void discardRandom(unsigned long long n, unsigned long long& this_random_count, std::mt19937& this_random_generator);
		unsigned long long getMaxRandomCount();
		void updateSingleSite(int i, int j, unsigned long long& this_random_count, std::mt19937& this_random_generator);

	private:

		// Lattice, dispersal, and time step parameters
		unsigned int seeds[5];
		int *species_abundance;
		std::mt19937 global_random_generator;
		std::string parameter_filename, outfile_base, outfile_dir;
		int id, num_species, lattice_size, num_steps, max_time_step, restart_time, min_persistence, num_restarts;
		unsigned long long random_count, max_random_count;
		double germination_probability, initial_occupancy;
		int **lattice, **next_lattice;
		double ***dispersal_lattice, ***next_dispersal_lattice;

		// Species specific parameters
		int *delta;
		double *species_occupancy, *juvenile_survival_probability, *adult_survival_probability;
		double *maximum_competition, *dispersal_probability, *dispersal_length, *intrinsic_fecundity;
	
		// Competition parameters
		std::string competition_filename, competition_type;
		double competition_lower_bound, competition_diag_lower_bound, competition_upper_bound, competition_diag_upper_bound; 
		double competition_mean, competition_diag_mean, competition_sdev, competition_diag_sdev;
		double competition_correlation, imbalance, fecundity_transitivity_type, growth_transitivity_type, fecundity_growth_relative_hierarchy;
		double fecundity_imbalance_mean, growth_imbalance_mean, fecundity_growth_correlation;
		double fecundity_relative_intransitivity, growth_relative_intransitivity;
		double **competition_fecundity, **competition_growth;
		double **fecundity_transitivity, **growth_transitivity; 
		double *fecundity_rank, *growth_rank;

		void initializeRandomSimulation();
		void initializeReplicateSimulation();
		void initializeRestartSimulation();
		void loadSeeds();
		void loadLattice();
		void loadDispersal();
		void loadCompetition();

		// Allocation and random seed
		std::mt19937& generateRandom();
		unsigned int getRandom();
		void allocSimulation();
		void initializeLattice();

		// Read input and set random parameters
		void checkInputFormat();
		void seedGenerator();
		void getParameter(int *value, std::string parameter_name, int essential);
		void getParameter(double *value, std::string parameter_name, int essential);
		void getParameter(std::string *value, std::string parameter_name, int essential);
		void getParameter(int *value_array, int n, std::string parameter_name, int essential);
		void getParameter(double *value_array, int n, std::string parameter_name, int essential);
		void initializeNormalRandomArray(double *vector, double *mean, double *sdev, int length);
		void setRandomParameter(double *intrinsic_fecundity, int num_species, std::string parameter_name, int type);
		std::string trimString(std::string str);
		std::string trimStringNoComment(std::string str);

		// Competition Initialization 
		void initializeUniformCompetition();
		void initializeTNormalCompetition();
		void initializeUniformCorrelatedCompetition();
		void initializeTNormalCorrelatedCompetition();
		void imbalanceCompetition();
		void setCompetitionTransitivity();
		void setGrowthCompetitionTransitivity( int * fecundity_hierarchy );

		// Competition properties
		void checkCorrelation();
		void getSpeciesAbundance();
		void getFecundityGrowthCorrelation();
		void getImbalanceMean();
		void getDiscreteTransitivity();
		void getDiscreteFecundityTransitivity();
		void getDiscreteGrowthTransitivity();
		void getContinuousTransitivity();

};

#endif
