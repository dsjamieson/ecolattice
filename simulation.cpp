// Percentage juvenile in initial conditions option
// optaional condiiton function for parameter file
// initialization defaults
// work out overfilling box problem

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <random>
#include <sys/stat.h>
#include "simulation.h"

Simulation::Simulation(std::string filename) {
	
	int i, j;

	//defaults
	delta = 1;
	competition_type = "Uniform";
	competition_upper_bound = 0;
	competition_lower_bound = -1;
	competition_correlation = 0;
	imbalance = 0.5;
	transitivity = 0.;
	parameter_filename = filename;

	// Read input
	getSeeds();
	getParameter(&box_size, "BoxSize", 1);
	getParameter(&num_species, "Species", 1);
	getParameter(&delta, "Delta", 0);
	getParameter(&max_time_step, "MaxTimeStep", 1);
	getParameter(&grid_occupancy, "GridOccupancy", 1);
	getParameter(&germination_probability, "GerminationProbability", 1);
	if(germination_probability > 1 || germination_probability < 0){
		fprintf(stdout, "GerminationProbability must be between 0 and 1\n");
		exit(0);
	}
	getParameter(&juvenile_survival_probability, "JuvenileSurvival", 1);
	if(juvenile_survival_probability > 1 || juvenile_survival_probability < 0){
		fprintf(stdout, "JuvenileSurvival must be between 0 and 1\n");
		exit(0);
	}
	getParameter(&adult_survival_probability, "AdultSurvival", 1);
	if(adult_survival_probability > 1 || adult_survival_probability < 0){
		fprintf(stdout, "AdultSurvival must be between 0 and 1\n");
		exit(0);
	}
	getParameter(&maximum_competition, "MaximumCompetition", 1);
	if(maximum_competition > 1 || maximum_competition < 0){
		fprintf(stdout, "MaximumCompetition must be between 0 and 1\n");
		exit(0);
	}
	getParameter(&dispersal_probability, "DispersalProbability", 1);
	if(dispersal_probability > 1 || dispersal_probability < 0){
		fprintf(stdout, "DispersalProbability must be between 0 and 1\n");
		exit(0);
	}

	species_occupancy = new double[num_species];
	dispersal_length = new double[num_species];
	intrinsic_fecundity = new double[num_species];
	getParameter(species_occupancy, num_species, "SpeciesOccupancy",  1);
	if(species_occupancy[0] == 1. && species_occupancy[num_species-1] == 1. ) {
		for(i=0; i<num_species; i++)
			species_occupancy[i] = species_occupancy[i]/(num_species);
	}	
	setRandomParameter(dispersal_length, num_species, "DispersalLength");
	setRandomParameter(intrinsic_fecundity, num_species, "Fecundity");

	getParameter(&competition_type, "CompetitionType", 0);
	getParameter(&competition_upper_bound, "CompetitionUpper", 0);
	getParameter(&competition_lower_bound, "CompetitionLower", 0);
	if(competition_type.compare("TNormal")==0)
		getParameter(&competition_sdev, "CompetitionSdev", 1);
	
	getParameter(&competition_correlation, "CompetitionCorr", 0);
	if(fabs(competition_correlation) > 1.) {
		fprintf(stderr, "Error, CompetitionCorr must be between -1 and 1\n");
		exit(0);
	}
	getParameter(&imbalance, "Imbalance", 0);
	getParameter(&transitivity, "Transitivity", 0);
	if( ( competition_correlation!=0 + imbalance!=0.5 + transitivity != 0 ) > 1  ) {
		fprintf(stderr, "Error, only one of CompetitionCorr, Imbalance, and Transitivty can be set\n");
		exit(0);
	}

	if(imbalance<0 || imbalance > 1) {
		fprintf(stderr, "Error, Imbalance must be between 0 and 1\n");
		exit(0);
	}

	if(transitivity*transitivity != 1. && transitivity != 0.  ) {
		fprintf(stderr, "Error, Transitivity must be -1, 0 or 1\n");
		exit(0);
	}
	getParameter(&outfile_base, "OutfileBase", 1);
	getParameter(&outfile_dir, "OutfileDir", 0);
	if(outfile_dir.size()!=0){
		struct stat buf;
		if(stat(outfile_dir.c_str(), &buf)==0)
			outfile_base = outfile_dir + "/" + outfile_base;
		else {
			fprintf(stdout, "Error, no directory %s\n", outfile_dir.c_str());
			exit(0);
		}
	}

	// Initializations
	if( competition_type.compare("Uniform") == 0  ) {
		if(transitivity!=0)
			initializeUniformTransitiveCompetition();
		else if(competition_correlation != 0) 
			initializeUniformCorrelatedCompetition();
		else
			initializeUniformCompetition();

	}
	else if( competition_type.compare("TNormal") == 0  ){
		if(transitivity!=0)
			initializeTNormalTransitiveCompetition();
		else if(competition_correlation != 0)
			initializeTNormalCorrelatedCompetition();
		else
			initializeTNormalCompetition();
	}
	else {
		fprintf(stderr, "Competition type %s not valid\n", competition_type.c_str());
		exit(0);
	}


	if(imbalance!=0.5)
		imbalanceCompetition();

	getImbalanceMean();
	getDiscreteTransitivity();
	getFecundityGrowthCorrelation();

	initializeBox();
	initializeNextBox();
	initializeDispersalBox();
	initializeNextDispersalBox();
	saveCompetition();
	saveProperties();

}

void Simulation::seedGenerator(int num_seeds) {

    std::seed_seq seq(seeds, seeds+num_seeds);
    std::vector<std::uint32_t> seeds(std::mt19937::state_size);
    seq.generate(seeds.begin(), seeds.end());
	global_random_generator.seed( seq );

	return;

}

int Simulation::getNewSeed(int i) {
	
	return i;

	}

void Simulation::initializeBox() {

	int i,j;

	// Initializes box with random distribution of the species
	std::bernoulli_distribution stage_dist(0.5);
	std::bernoulli_distribution occupy_dist(grid_occupancy);
	std::discrete_distribution<int> species_dist(species_occupancy, species_occupancy+num_species);


	box = new int *[box_size];
	for( i = 0 ; i < box_size ; i++ ){
		box[i] = new int[box_size];
		for(j=0;j<box_size;j++) {
			box[i][j] = occupy_dist(global_random_generator);
			if(stage_dist(global_random_generator) && box[i][j] )
				box[i][j] *= -1;
			if(box[i][j] != 0)
				box[i][j] *=1+species_dist(global_random_generator);
		}
	}

	return;

}

void Simulation::initializeNextBox() {

	int i,j;

	// Initializes box with random distribution of the species
	//
	next_box = new int *[box_size];
	for( i = 0 ; i < box_size ; i++ ){
		next_box[i] = new int[box_size];
		for(j=0;j<box_size;j++) {
			next_box[i][j] = 0;
		}
	}

	return;

}



void Simulation::initializeDispersalBox() {

	int i,j,k;

	// Initializes box with random distribution of the species
	dispersal_box = new double **[box_size];
	for( i = 0 ; i < box_size ; i++ ){
		dispersal_box[i] = new double *[box_size];
		for(j=0;j<box_size;j++) {
			dispersal_box[i][j] = new double[num_species];
			for(k=0;k<num_species;k++)
				dispersal_box[i][j][k] = 0;

		}
	}

	return;

}


void Simulation::initializeNextDispersalBox() {

	int i,j,k;

	// Initializes box with random distribution of the species
	next_dispersal_box = new double **[box_size];
	for( i = 0 ; i < box_size ; i++ ){
		next_dispersal_box[i] = new double *[box_size];
		for(j=0;j<box_size;j++) {
			next_dispersal_box[i][j] = new double[num_species];
			for(k=0;k<num_species;k++)
				next_dispersal_box[i][j][k] = 0;

		}
	}

	return;

}




void Simulation::updateSingleSite(int i, int j, std::mt19937& local_random_generator) {


	int k, l, kp, lp;

	int this_species = abs(box[i][j]);
	int this_stage;
	double growth_probability = 0;
	int total_abundance = 0;
	double this_fecundity = 0.;

	int neighborhood[ (2*delta+1)*(2*delta+1) ]; 
	int neighborhood_abundance[num_species];
	for(k=0;k<(2*delta+1)*(2*delta+1);k++)
		neighborhood[k]=0;
	for(k=0;k<num_species;k++)
		neighborhood_abundance[k]=0;

	double **distance_probability;
	double distance_probability_sum = 0.;

	if(this_species==0) {

		int total_seeds = 0;
	
		for(k=0;k<num_species;k++)
			total_seeds += dispersal_box[i][j][k];

		std::bernoulli_distribution germdist(1. - pow(1.- germination_probability, total_seeds ));
		
		if(germdist(local_random_generator)) {
			std::discrete_distribution<int> species_dist(dispersal_box[i][j], dispersal_box[i][j]+num_species);
			next_box[i][j] = - species_dist(local_random_generator);
		}

	}
	else {

		this_stage = box[i][j]/abs(box[i][j]);

		double survival_probability;
	
		if(this_stage == 1)
			survival_probability = adult_survival_probability ;
		else
			survival_probability = juvenile_survival_probability;

		std::bernoulli_distribution survival_dist(survival_probability);

		if(survival_dist(local_random_generator)) {

			next_box[i][j] = box[i][j];	

			for(k=0; k<2*delta+1; k++) {
				for(l=0;l<2*delta+1;l++) {
					if(k!=delta || l!=delta ){

						kp = k;
						lp = l;

						if( i - delta + kp < 0)
							kp += box_size - i;
						else if( i - delta + kp > box_size-1)
							kp += -box_size;	
						if( j - delta + lp < 0)
							lp += box_size - j;
						else if( j - delta + lp > box_size-1)
							lp += -box_size;

						neighborhood[l+(2*delta+1)*k] = box[i-delta+kp][j-delta+lp];
						if( neighborhood[l+(2*delta+1)*k] != 0  ) {
							neighborhood_abundance[ abs( neighborhood[l+(2*delta+1)*k] )  - 1  ]++;
							total_abundance++;
						}
						
					}

				}

			}

			if( this_stage < 0  ) {
		
			for(k=0; k<num_species; k++)
				growth_probability += competition_growth[this_species-1][k]*( (double) neighborhood_abundance[k] )/( (double) total_abundance );
			growth_probability = maximum_competition*exp(growth_probability);
			if(growth_probability > maximum_competition)
				growth_probability = maximum_competition;

			std::bernoulli_distribution growth_dist(growth_probability);

			if(growth_dist(local_random_generator))
				next_box[i][j] *= -1;
			}
			else {



				for(k=0; k<num_species; k++)
					this_fecundity += competition_fecundity[this_species-1][k]*( (double) neighborhood_abundance[k])/( (double) total_abundance );
				this_fecundity = intrinsic_fecundity[this_species-1]*exp(this_fecundity);



				distance_probability = new double*[box_size];
				for(k=0; k<box_size;k++) {
					distance_probability[k] = new double[box_size];
					for(l=0;l<box_size;l++) {
						distance_probability[k][l] = pow(std::min(abs(i-k), box_size -  abs(i-k)),2);
						distance_probability[k][l] += pow(std::min(abs(j-l), box_size -  abs(j-l)),2);
						distance_probability[k][l] = sqrt(distance_probability[k][l]);
						distance_probability[k][l] = exp(log(dispersal_probability)/dispersal_length[this_species-1]*distance_probability[k][l]);
						distance_probability_sum += distance_probability[k][l];

					}
				}

				for(k=0; k<box_size;k++) {
					for(l=0;l<box_size;l++) {
						next_dispersal_box[k][l][this_species-1] = this_fecundity*distance_probability[k][l]/distance_probability_sum;
					}
				}
			}
	

		}
		else {
			next_box[i][j] = 0;
		}

	}
	
	return;

}

void Simulation::nextToThis() {
	
	int i, j, k;

	for(i=0;i<box_size;i++) {
		for(j=0;j<box_size;j++) {

			box[i][j] = next_box[i][j];
			next_box[i][j] = 0;
	
			for(k=0;k<num_species;k++) {
				dispersal_box[i][j][k] = next_dispersal_box[i][j][k];
				next_dispersal_box[i][j][k];
			}
		}
	}

	return;

}

void Simulation::saveBox(int time_step) {

	int i, j;

	

	std::ofstream boxfile;
	boxfile.open(outfile_base+"_"+std::to_string(time_step)+".csv", std::ios::out | std::ios::trunc);

	for(i=0;i<box_size;i++) {
		for(j=0;j<box_size;j++) {

			if(box[i][j] < 0)
				boxfile << box[i][j];
			else
				boxfile << " " << box[i][j];

			if(j!=box_size-1)
				boxfile << ", ";

			}
			boxfile << std::endl;
		}

	boxfile.close();

	return;

}

int Simulation::getBoxSize() {
	return box_size;
}

int Simulation::getMaxTimeStep() {
	return max_time_step;
}

void Simulation::saveCompetition() {

	int i, j;
	
	std::ofstream competition_file;
	competition_file.open(outfile_base+"_competition.csv", std::ios::out | std::ios::trunc);


	competition_file << "# Intrinsic fecundity:" << std::endl;
	for(i=0;i<num_species;i++) {
		competition_file << " " << intrinsic_fecundity[i];
	}
	competition_file << std::endl;

	competition_file << "# Dispersal lengths:" << std::endl;
	for(i=0;i<num_species;i++) {
		competition_file << " " << intrinsic_fecundity[i];
	}
	competition_file << std::endl;

	competition_file << "# Fecundity competition:" << std::endl;
	for(i=0;i<num_species;i++) {
		for(j=0;j<num_species;j++) {

			if(competition_fecundity[i][j] < 0)
				competition_file << competition_fecundity[i][j];
			else
				competition_file << " " << competition_fecundity[i][j];

			if(j!=num_species-1)
				competition_file << ", ";

			}
			competition_file << std::endl;
		}

	competition_file << "# Growth competition:" << std::endl;
	for(i=0;i<num_species;i++) {
		for(j=0;j<num_species;j++) {

			if(competition_growth[i][j] < 0)
				competition_file << competition_growth[i][j];
			else
				competition_file << " " << competition_growth[i][j];

			if(j!=num_species-1)
				competition_file << ", ";

			}
			competition_file << std::endl;
		}

	competition_file.close();

	return;

}

void Simulation::saveProperties() {

	std::ofstream properties_file;
	properties_file.open(outfile_base+"_properties.csv", std::ios::out | std::ios::trunc);

	properties_file <<  "# Fecundity imbalance mean:" << std::endl  << fecundity_imbalance_mean <<  std::endl;
	properties_file <<  "# Growth imbalance mean:" << std::endl  << growth_imbalance_mean <<  std::endl;
	properties_file <<  "# Fecundity relative intransitivity:" << std::endl  << fecundity_relative_intransitivity <<  std::endl;
	properties_file <<  "# Growth relative intransitivity:" << std::endl  << growth_relative_intransitivity <<  std::endl;
	properties_file <<  "# Fecundity-growth cross correlation:" << std::endl  << fecundity_growth_correlation <<  std::endl;

	properties_file.close();

	return;

}
