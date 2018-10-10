#ifndef SIMULATION_H
#define SIMULATION_H

extern std::mt19937 global_random_generator;

class Simulation {

	public:
		Simulation(std::string filename);
		int getNewSeed(int i);
		int getBoxSize();
		int getMaxTimeStep();
		void updateSingleSite(int i, int j, std::mt19937& local_random_generator);
		void saveBox(int time_step);
		void nextToThis();


	private:
		int num_species, box_size, num_steps, delta, seed, max_time_step;
		double grid_occupancy, imbalance, transitivity, fecundity_imbalance_mean, growth_imbalance_mean;
		double fecundity_relative_intransitivity, growth_relative_intransitivity, fecundity_growth_correlation;
		double germination_probability, juvenile_survival_probability, adult_survival_probability;
		double maximum_competition, dispersal_probability;
		std::string competition_type, outfile_base, outfile_dir;
		double *intrinsic_fecundity, *intrinsic_fecundity_sdev, *dispersal_length, *species_occupancy;
		double **competition_fecundity, **competition_growth;
		double **fecundity_transitivity, **growth_transitivity ; 
		double *fecundity_row_sum, *growth_row_sum;
		double competition_lower_bound, competition_upper_bound, competition_sdev, competition_correlation;
		int **box, **next_box;
		double ***dispersal_box, ***next_dispersal_box;

		void seedGenerator();
		std::string trimstr(std::string str);

		void getParameter(int *value, std::string parameter_name, std::string filename, int essential);
		void getParameter(double *value, std::string parameter_name, std::string filename, int essential);
		void getParameter(double *value_array, int n, std::string parameter_name, std::string filename, int essential);
		void getParameter(std::string *value, std::string parameter_name, std::string filename, int essential);
		void initializeNormalRandomArray(double *vector, double *mean, double *sdev, int length);
		void setRandomParameter(double *intrinsic_fecundity, int num_species, std::string parameter_name, std::string parameter_filename);

		void initializeUniformCompetition();
		void initializeTNormalCompetition();
		void initializeUniformCorrelatedCompetition();
		void initializeTNormalCorrelatedCompetition();
    	void initializeUniformTransitiveCompetition();
   		void initializeTNormalTransitiveCompetition();
		void imbalanceCompetition();
		void checkCorrelation();

		void getFecundityGrowthCorrelation();
		void getImbalanceMean();
		void getDiscreteTransitivity();
		void getDiscreteFecundityTransitivity();
		void getContinuousTransitivity();

		void initializeBox();
		void initializeNextBox();
		void initializeDispersalBox();
		void initializeNextDispersalBox();

		void saveCompetition();
		void saveProperties();


};

#endif
