#ifndef SIMULATION_H
#define SIMULATION_H

#include <string>
#include <random>
#include <fstream>
#include <sstream>
#include <cmath>
#include <mpi.h>
#include <sys/stat.h>
#include "simulation.h"

class Simulation {

	public:
		Simulation(std::string filename, int p_id);
		int getBoxSize();
		int getSpecies();
		int getMaxTimeStep();
		int getSite(int i, int j);
		int getNextSite(int i, int j);
		double getDispersal(int i, int j, int s);
		double getNextDispersal(int i, int j, int s);
		int getRestartTime();


		void setSite(int i, int j, int s);
		void addSite(int i, int j, int s);
		void setDispersal(int i, int j, int s, double t);
		void addDispersal(int i, int j, int s, double t);
		void resetBox();
		void resetNextBox();

		void updateSingleSite(int i, int j);
		void discardRandom(unsigned long long n);
		unsigned long long getRandomCount();
		void saveBox(int time_step);
		void saveCompetition();
		void saveProperties();
		void saveDispersal(int time_step);
		void nextToThis();



	private:

		// Box, dispersal, and time step parameters
		int *seeds;
		std::mt19937 global_random_generator;
		std::string parameter_filename, outfile_base, outfile_dir;
		int id, num_species, box_size, num_steps, max_time_step, restart_time;
		unsigned long long random_count, max_random_count;
		double germination_probability, initial_occupancy;
		int **box, **next_box;
		double ***dispersal_box, ***next_dispersal_box;

		// Species specific parameters
		int *delta;
		double *species_occupancy, *juvenile_survival_probability, *adult_survival_probability;
		double *maximum_competition, *dispersal_probability, *dispersal_length, *intrinsic_fecundity;
	
		// Competition parameters
		std::string competition_filename, competition_type;
		double competition_lower_bound, competition_upper_bound, competition_mean, competition_sdev;
		double competition_correlation, imbalance, fecundity_transitivity_type, growth_transitivity_type, fecundity_growth_relative_hierarchy;
		double fecundity_imbalance_mean, growth_imbalance_mean, fecundity_growth_correlation;
		double fecundity_relative_intransitivity, growth_relative_intransitivity;
		double **competition_fecundity, **competition_growth;
		double **fecundity_transitivity, **growth_transitivity; 
		double *fecundity_row_sum, *growth_row_sum;

		void initializeRandomSimulation();
		void initializeRedoSimulation();
		void initializeRestartSimulation();

		// Allocation and random seed
		std::mt19937& generateRandom();
		int getRandom();
		void allocSim();
		void initializeBox();

		// Read input and set random parameters
		void checkInputFormat();
		void getSeeds();
		void seedGenerator(int num_seeds);
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
		void getFecundityGrowthCorrelation();
		void getImbalanceMean();
		void getDiscreteTransitivity();
		void getDiscreteFecundityTransitivity();
		void getDiscreteGrowthTransitivity();
		void getContinuousTransitivity();

		void loadBox();
		void loadDispersal();
		void loadCompetition();

};

#endif
