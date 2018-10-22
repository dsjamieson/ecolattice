
	 /**********************************************************
	 *
	 *			D.S Jamieson and N.L Kinlock, 2018			
	 *
	 *	The Simulation class constructor is defined here, 
	 *	requiring an imput file containing the simulation 
	 *	parameters. Public methods for getting, setting and
	 *	adding to values of private arrays and variable are 
	 *	here, as well as the methods for saving output.
	 *
	 ***********************************************************/

#include "simulation.h"

Simulation::Simulation(std::string filename, int p_id) {
	
	id = p_id;
	random_count = 0;

	//defaults
	restart_time = 0;
	delta = 1;
	competition_type = "Uniform";
	competition_upper_bound = 0;
	competition_lower_bound = -1;
	competition_correlation = 0;
	imbalance = 0.5;
	fecundity_transitivity_type = 0.;
	growth_transitivity_type = 0.;
	fecundity_growth_relative_hierarchy = 0.;
	parameter_filename = filename;

	// Checks for any obvious format issues with input parameter file 
	checkInputFormat();

	// Checks if this is continuing a previous simulation from time step retart_time
	getParameter(&restart_time, "RestartTime", 0);
	if( restart_time < 0) {
		if(id==0)
			fprintf(stderr, "Error, RestartTime must be positive\n");
		MPI_Finalize();
		exit(0);
	}

	if(restart_time != 0) {
		getParameter(&competition_filename, "CompetitionFile", 1);
	}
	else {
		getParameter(&competition_filename, "CompetitionFile", 0);
	}

	// Simulation parameters
	getParameter(&box_size, "BoxSize", 1);
	if( box_size < 0) {
		if(id==0)
			fprintf(stderr, "Error, BoxSize must be positive\n");
		MPI_Finalize();
		exit(0);
	}

	getParameter(&num_species, "Species", 1);
	getParameter(&delta, "Delta", 0);
	if(delta < 0 || delta > box_size - 1) {
		if(id==0)
			fprintf(stderr, "Error, Delta must be greater than zero and less than BoxSize\n");
		MPI_Finalize();
		exit(0);
	}
	getParameter(&max_time_step, "MaxTimeStep", 1);
	getParameter(&initial_occupancy, "InitialOccupancy", 1);
	if(initial_occupancy > 1 || initial_occupancy < 0){
		if(id==0)
			fprintf(stderr, "Error, InitialOccupancy must be between 0 and 1\n");
		MPI_Finalize();
		exit(0);
	}
	getParameter(&germination_probability, "GerminationProbability", 1);
	if(initial_occupancy > 1 || initial_occupancy < 0){
		if(id==0)
			fprintf(stderr, "Error, GerminationProbability must be between 0 and 1\n");
		MPI_Finalize();
		exit(0);
	}
	getParameter(&outfile_base, "OutfileBase", 1);
	getParameter(&outfile_dir, "OutfileDir", 0);
	if(outfile_dir.size()!=0){
		struct stat buf;
		if(stat(outfile_dir.c_str(), &buf)==0)
			outfile_base = outfile_dir + "/" + outfile_base;
		else {
			if(id==0)
				fprintf(stderr, "Error, no directory %s\n", outfile_dir.c_str());
			MPI_Finalize();
			exit(0);
		}
	}


	// Allocate simulation arrays and seed random generator
	max_random_count = (int64_t) 1000.*(4.*num_species*num_species + 5.*box_size*box_size);
	allocSim();
	getSeeds();

	
	if( restart_time != 0) {
		// Continue previous simulation
		initializeRestartSimulation();
	}
	else if( competition_filename.size() != 0  ) {
		// Use specific competition parameters
		initializeRedoSimulation();
	}
	else {
		// Start a new simulation
		initializeRandomSimulation();
	}


}


void Simulation::initializeRandomSimulation() {

	initializeBox();

	// Set species specific parameters
	setRandomProbability(species_occupancy, num_species, "SpeciesOccupancy", 3);
	setRandomProbability(juvenile_survival_probability, num_species, "JuvenileSurvival", 2);
	setRandomProbability(adult_survival_probability, num_species, "AdultSurvival", 2);
	setRandomProbability(maximum_competition, num_species, "MaximumCompetition", 2);
	setRandomProbability(dispersal_probability,  num_species,"DispersalProbability", 2);
	setRandomParameter(dispersal_length, num_species, "DispersalLength");
	setRandomParameter(intrinsic_fecundity, num_species, "Fecundity");

	// Set competition parameters
	getParameter(&competition_lower_bound, "CompetitionLower", 0);
	if(fabs(competition_lower_bound) > 1.) {
		if(id==0)
			fprintf(stderr, "Error, CompetitionLower must be between -1 and 1\n");
		MPI_Finalize();
		exit(0);
	}
	getParameter(&competition_upper_bound, "CompetitionUpper", 0);
	if(fabs(competition_upper_bound) > 1.) {
		if(id==0)
			fprintf(stderr, "Error, CompetitionUpper must be between -1 and 1\n");
		MPI_Finalize();
		exit(0);
	}
	competition_mean = (competition_lower_bound+competition_upper_bound)/2.;

	getParameter(&competition_type, "CompetitionType", 0);
	if(competition_type.compare("TNormal")==0) {
		getParameter(&competition_mean, "CompetitionMean", 0);
		if(competition_mean < competition_lower_bound || competition_mean > competition_upper_bound) {
			if(id==0)
				fprintf(stderr, "Error, CompetitionMean must be between CompetitionLower and CompetitionUpper\n");
			MPI_Finalize();
			exit(0);
		}
		getParameter(&competition_sdev, "CompetitionSdev", 1);
	}

	getParameter(&competition_correlation, "CompetitionCorr", 0);
	if(fabs(competition_correlation) > 1.) {
		if(id==0)
			fprintf(stderr, "Error, CompetitionCorr must be between -1 and 1\n");
		MPI_Finalize();
		exit(0);
	}

	getParameter(&imbalance, "Imbalance", 0);
	if( imbalance < 0 || imbalance > 1) {
		if(id==0)
			fprintf(stderr, "Error, Imbalance must be between 0 and 1\n");
		MPI_Finalize();
		exit(0);
	}

	getParameter(&fecundity_transitivity_type, "FecundityTransitivity", 0);
	if( fabs(fecundity_transitivity_type) != 1. && fecundity_transitivity_type != 0.) {
		if(id==0)
			fprintf(stderr, "Error, current implementation only allows for FecunidtyTransitivity 0 (random), maximum (1) or minimum (-1)\n");
		MPI_Finalize();
		exit(0);
	}
	getParameter(&growth_transitivity_type, "GrowthTransitivity", 0);
	if( fabs(growth_transitivity_type) != 1. && growth_transitivity_type != 0.) {
		if(id==0)
			fprintf(stderr, "Error, current implementation only allows for GrowthTransitivity 0 (random), maximum (1) or minimum (-1)\n");
		MPI_Finalize();
		exit(0);
	}
	getParameter(&fecundity_growth_relative_hierarchy, "RelativeHierarchy", 0);
	if( fabs(fecundity_growth_relative_hierarchy) != 1. && fecundity_growth_relative_hierarchy != 0.) {
		if(id==0)
			fprintf(stderr, "Error, RelativeHierarchy must be +/-1 (equal/univerted), or 0 (unrelated)\n");
		MPI_Finalize();
		exit(0);
	}
	if( fecundity_growth_relative_hierarchy != 0 && ( fecundity_transitivity_type == 0 ||  growth_transitivity_type == 0  )  ) {
		if(id==0)
			fprintf(stderr, "Error, if RelativeHierarchy is not zero, neither FecundityTransitive nor GrowthTransitivity can be zero\n");
		MPI_Finalize();
		exit(0);
	}
	if( ( competition_correlation!=0 ) + ( imbalance!=0.5 ) + ( ( fabs(fecundity_transitivity_type) + fabs(growth_transitivity_type) ) != 0 ) > 1  ) {
		if(id==0)
			fprintf(stderr, "Error, only one of CompetitionCorr, Imbalance, and (Fecundity/Growth)Transitivty can be set\n");
		MPI_Finalize();
		exit(0);
	}

	// Initialization competition
	if( competition_type.compare("Uniform") == 0  ) {
		if(competition_correlation != 0) 
			initializeUniformCorrelatedCompetition();
		else
			initializeUniformCompetition();

	}
	else if( competition_type.compare("TNormal") == 0  ){
		if(competition_correlation != 0)
			initializeTNormalCorrelatedCompetition();
		else
			initializeTNormalCompetition();
	}
	else {
		if(id==0)
			fprintf(stderr, "Competition type must be either Uniform or TNormal\n", competition_type.c_str());
		MPI_Finalize();
		exit(0);
	}
	if(imbalance!=0.5)
		imbalanceCompetition();

	if(fecundity_transitivity_type!=0 || growth_transitivity_type!=0)
		setCompetitionTransitivity();

	// Calculate competition properties
	getImbalanceMean();
	getDiscreteTransitivity();
	getFecundityGrowthCorrelation();

	if(random_count > max_random_count) {
		if(id==0)
			fprintf(stderr, "Error, too many random numbers used to generate initial conditions, for TNormal distribution or correlation may be too sever\n");
		MPI_Finalize();
		exit(0);
	}

	global_random_generator.discard(max_random_count - random_count);
	random_count = 0;

	return;

}

void Simulation::initializeRedoSimulation() {

	initializeBox();
	loadCompetition();

	if(random_count > max_random_count) {
		if(id==0)
			fprintf(stderr, "Error, too many random numbers used to generate initial conditions, for TNormal distribution or correlation may be too sever\n");
		MPI_Finalize();
		exit(0);
	}

	global_random_generator.discard(max_random_count - random_count);
	random_count = 0;

	return;

}

void Simulation::initializeRestartSimulation() {

	loadBox();
	loadDispersal();
	loadCompetition();

	if(random_count > max_random_count) {
		if(id==0)
			fprintf(stderr, "Error, too many random numbers used to generate initial conditions, for TNormal distribution or correlation may be too sever\n");
		MPI_Finalize();
		exit(0);
	}

	global_random_generator.discard(max_random_count + 4*box_size*box_size*restart_time - random_count);
	random_count = 0;

	return;

}

void Simulation::allocSim() {

	int i, j, k;

	// Allocate simulation arrays
	juvenile_survival_probability = new double[num_species];
	adult_survival_probability = new double[num_species];
	maximum_competition = new double[num_species];
	dispersal_probability = new double[num_species];
	if( !juvenile_survival_probability || !adult_survival_probability || !maximum_competition) {
		fprintf(stderr, "Error, unable to allocate memory for species parameters\n");
		MPI_Finalize();
		exit(-1);
	}
	fecundity_row_sum = new double[num_species];
	growth_row_sum = new double[num_species];
	fecundity_transitivity = new double *[num_species];
	growth_transitivity = new double *[num_species];
	if(!fecundity_row_sum || !growth_row_sum || !fecundity_transitivity || !growth_transitivity) {
		fprintf(stderr, "Error, unable to allocate memory for transitivity arrays 1\n");
		MPI_Finalize();
		exit(-1);
	}
	species_occupancy = new double[num_species];
	dispersal_length = new double[num_species];
	intrinsic_fecundity = new double[num_species];
	if(!species_occupancy  || !dispersal_length|| !intrinsic_fecundity ) {
		fprintf(stderr, "Error, unable to allocate memory for species parameters\n");
		MPI_Finalize();
		exit(-1);
	}
	competition_fecundity = new double *[num_species];
	competition_growth = new double *[num_species];
	if(! competition_fecundity || ! competition_growth) {
		fprintf(stderr, "Error, unable to allocate memory for competition matrices\n");
		MPI_Finalize();
		exit(-1);
	}
	for( i = 0 ; i < num_species ; i++ ){

		juvenile_survival_probability[i] = 0.;
		adult_survival_probability[i] = 0.;
		maximum_competition[i] = 0.;
		dispersal_probability[i] = 0.;
		species_occupancy[i] = 0.;
		dispersal_length[i] = 0.;
		intrinsic_fecundity[i] = 0.;
		fecundity_row_sum[i] = 0.;
		growth_row_sum[i] = 0.;

		competition_fecundity[i] = new double[num_species];
		competition_growth[i] = new double[num_species];
		fecundity_transitivity[i] = new double[num_species];
		growth_transitivity[i] = new double[num_species];
		if(!fecundity_transitivity[i] || !growth_transitivity[i]  ) {
			fprintf(stderr, "Error, unable to allocate memory for transitivity arrays\n");
			MPI_Finalize();
			exit(-1);
		}
		if(!competition_fecundity[i] || !competition_growth[i]) {
			fprintf(stderr, "Error, unable to allocate memory for competition matrices\n");
			MPI_Finalize();
			exit(-1);
		}
		for(j=0;j<num_species;j++) {
			fecundity_transitivity[i][j] = 0.;
			growth_transitivity[i][j] = 0.;
			competition_fecundity[i][j] = 0.;
			competition_growth[i][j] = 0.;
		}

	}

	box = new int *[box_size];
	next_box = new int *[box_size];
	dispersal_box = new double **[box_size];
	next_dispersal_box = new double **[box_size];
	if(!box  || !next_box || !dispersal_box || !next_dispersal_box) {
		fprintf(stderr, "Error, unable to allocate memory for box and dispersal\n");
		MPI_Finalize();
		exit(-1);
	}
	for( i = 0 ; i < box_size ; i++ ){

		box[i] = new int[box_size];
		next_box[i] = new int[box_size];
		dispersal_box[i] = new double *[box_size];
		next_dispersal_box[i] = new double *[box_size];
		if(!box[i]  || !next_box[i] || !dispersal_box[i] || !next_dispersal_box[i]) {
			fprintf(stderr, "Error, unable to allocate memory for box and dispersal\n");
			MPI_Finalize();
			exit(-1);
		}

		for(j=0;j<box_size;j++) {

			box[i][j] = 0;
			next_box[i][j] = 0;

			dispersal_box[i][j] = new double[num_species];
			next_dispersal_box[i][j] = new double[num_species];
			if( !dispersal_box[i][j] || !next_dispersal_box[i][j]) {
				fprintf(stderr, "Error, unable to allocate memory for dispersal\n");
				MPI_Finalize();
				exit(-1);
			}


			for(k=0;k<num_species;k++) {
				dispersal_box[i][j][k] = 0;
				next_dispersal_box[i][j][k] = 0;
			}

		}
	}

	return;

}


void Simulation::initializeBox() {

	int i,j;

	// Initializes box with random distribution of the species
	std::bernoulli_distribution stage_dist(0.5);
	std::bernoulli_distribution occupy_dist(initial_occupancy);
	std::discrete_distribution<int> species_dist(species_occupancy, species_occupancy+num_species);

	for( i = 0 ; i < box_size ; i++ ){
		for(j=0;j<box_size;j++) {

			box[i][j] = occupy_dist(generateRandom());

			if( stage_dist(generateRandom()) == 0 && box[i][j] != 0 )
				box[i][j] *= -1;

			if(box[i][j] != 0) {
				box[i][j] *= 1+species_dist(generateRandom());
			}

		}
	}

	return;

}


void Simulation::resetBox() {
	
	int i, j, k;

	for(i=0;i<box_size;i++) {
		for(j=0;j<box_size;j++) {
			box[i][j] = 0;
			for(k=0;k<num_species;k++) {
				dispersal_box[i][j][k] = 0.;
			}
		}
	}

	return;

}


void Simulation::resetNextBox() {
	
	int i, j, k;

	for(i=0;i<box_size;i++) {
		for(j=0;j<box_size;j++) {
			next_box[i][j] = 0;
			for(k=0;k<num_species;k++) {
				next_dispersal_box[i][j][k] = 0.;
			}
		}
	}

	return;

}

std::mt19937& Simulation::generateRandom() {

	random_count+=2;
	return global_random_generator;

}


void Simulation::updateSingleSite(int i, int j) {

	int k, l, kp, lp;
	unsigned long long start_random_count = random_count;

	int this_species = 0;
	int this_stage = 0;
	double this_survival_probability = 0;
	double this_intrinsic_fecundity = 0.;
	double this_dispersal_probability = 0.;
	double this_dispersal_length = 0.;
	double this_maximum_competition = 0.;

	double total_seeds = 0;
	int total_abundance = 0;

	double growth_probability = 0;
	double this_fecundity = 0.;
	double distance_probability_sum = 0.;

	int *neighborhood;
	int *neighborhood_abundance;
	double **distance_probability;

	this_species = abs(box[i][j]);

	if(this_species==0) {

		total_seeds = 0;
	
		for(k=0;k<num_species;k++)
			total_seeds += dispersal_box[i][j][k];

		std::bernoulli_distribution germ_dist(1. - pow(1.- germination_probability, total_seeds ));

		if(germ_dist(generateRandom())) {
			std::discrete_distribution<int> species_dist(dispersal_box[i][j], dispersal_box[i][j]+num_species);
			next_box[i][j] = -abs( species_dist(generateRandom()) + 1 );
		}
		else {
			next_box[i][j] = 0;
		}

	}
	else {

		this_stage = box[i][j]/abs(box[i][j]);
	
		if(this_stage > 0)
			this_survival_probability = adult_survival_probability[this_species-1] ;
		else
			this_survival_probability = juvenile_survival_probability[this_species-1];

		this_intrinsic_fecundity = intrinsic_fecundity[this_species-1];
		this_dispersal_probability = dispersal_probability[this_species-1];
		this_dispersal_length = dispersal_length[this_species-1];
		this_maximum_competition = maximum_competition[this_species-1];

		std::bernoulli_distribution survival_dist(this_survival_probability);

		if(survival_dist(generateRandom())) {

			next_box[i][j] = box[i][j];

			neighborhood = new int[ (2*delta+1)*(2*delta+1) ]; 
			if(!neighborhood){
					fprintf(stderr, "Error, neighborhood memory allocation failed for site %d %d\n", i, j);	
					MPI_Finalize();
					exit(-1);
			}
			for(k=0;k<(2*delta+1)*(2*delta+1);k++)
				neighborhood[k]=0;

			neighborhood_abundance = new int[num_species];
			if(!neighborhood_abundance){
					fprintf(stderr, "Error, neighborhood abundance memory allocation failed for site %d %d\n", i, j);	
					MPI_Finalize();
					exit(-1);
			}
			for(k=0;k<num_species;k++)
				neighborhood_abundance[k]=0;	

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
		
			delete[] neighborhood;

			if( this_stage < 0  ) {

				for(k=0; k<num_species; k++)
					growth_probability += competition_growth[this_species-1][k]*( (double) neighborhood_abundance[k] )/( (double) total_abundance );
				growth_probability = this_maximum_competition*exp(growth_probability);
				if(growth_probability > this_maximum_competition)
					growth_probability = this_maximum_competition;

				std::bernoulli_distribution stage_dist(growth_probability);

				if(stage_dist(generateRandom()))
					next_box[i][j] = abs(next_box[i][j]);

			}
			else {

				if(total_abundance != 0) {
					for(k=0; k<num_species; k++)
						this_fecundity += competition_fecundity[this_species-1][k]*( (double) neighborhood_abundance[k])/( (double) total_abundance );
					this_fecundity = this_intrinsic_fecundity*exp(this_fecundity);
				}
				else {
					this_fecundity = this_intrinsic_fecundity;
				}


				distance_probability = new double*[box_size];
				if(!distance_probability) {
					fprintf(stderr, "Error, distance probability memory allocation failed for site %d %d\n", i, j);	
					MPI_Finalize();
					exit(-1);
				}
				for(k=0; k<box_size;k++) {
					distance_probability[k] = new double[box_size];
					if(!distance_probability[k]) {
						fprintf(stderr, "Error, distance probability memory allocation failed for site %d %d\n", i, j);	
						MPI_Finalize();
						exit(-1);
					}
					for(l=0;l<box_size;l++) {
						distance_probability[k][l] = pow(std::min(abs(i-k), box_size -  abs(i-k)),2);
						distance_probability[k][l] += pow(std::min(abs(j-l), box_size -  abs(j-l)),2);
						distance_probability[k][l] = sqrt(distance_probability[k][l]);
						distance_probability[k][l] = exp(log(this_dispersal_probability)/dispersal_length[this_species-1]*distance_probability[k][l]);
						distance_probability_sum += distance_probability[k][l];
					}
				}

				

				for(k=0; k<box_size;k++) {
					for(l=0;l<box_size;l++) {
						next_dispersal_box[k][l][this_species-1] += this_fecundity*distance_probability[k][l]/distance_probability_sum;
					}

					delete[] distance_probability[k];

				}

				delete[] distance_probability;

			}
			
			delete[] neighborhood_abundance;

		}
		else {

			next_box[i][j] = 0;
		}

	}

	discardRandom( ( (unsigned long long)  4 - ( random_count - start_random_count)  ) );

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
				next_dispersal_box[i][j][k] = 0;
			}
		}
	}

	return;

}


void Simulation::saveBox(int time_step) {

	int i, j;	

	std::ofstream box_file;
	box_file.open(outfile_base+"_"+std::to_string(time_step)+".csv", std::ios::out | std::ios::trunc);

	if(!box_file.is_open()) {
			if(id==0)
				fprintf(stderr, "Error, could not open time step %d box file file for for loading\n", time_step);
			MPI_Finalize();
			exit(0);
	}


	for(i=0;i<box_size;i++) {
		for(j=0;j<box_size;j++) {

			if(box[i][j] < 0)
				box_file << box[i][j];
			else
				box_file << " " << box[i][j];

			if(j!=box_size-1)
				box_file << ", ";

			}
			box_file << std::endl;
		}

	box_file.close();

	return;

}


void Simulation::saveCompetition() {

	int i, j;

	std::ofstream competition_file;
	competition_file.open(outfile_base+"_competition.csv", std::ios::out | std::ios::trunc);

	if(!competition_file.is_open()) {
			if(id==0)
				fprintf(stderr, "Error, could not open competition file for for saving\n");
			MPI_Finalize();
			exit(0);
	}


	competition_file << "# Species Occupancy:" << std::endl;
	for(i=0;i<num_species;i++) {
		competition_file << " " << species_occupancy[i];
	}
	competition_file << std::endl;

	competition_file << "# Juvenile Survival:" << std::endl;
	for(i=0;i<num_species;i++) {
		competition_file << " " << juvenile_survival_probability[i];
	}
	competition_file << std::endl;

	competition_file << "# Adult Survival:" << std::endl;
	for(i=0;i<num_species;i++) {
		competition_file << " " << adult_survival_probability[i];
	}
	competition_file << std::endl;

	competition_file << "# Maximum Competition:" << std::endl;
	for(i=0;i<num_species;i++) {
		competition_file << " " << maximum_competition[i];
	}
	competition_file << std::endl;

	competition_file << "# Dispersal probability:" << std::endl;
	for(i=0;i<num_species;i++) {
		competition_file << " " << dispersal_probability[i];
	}
	competition_file << std::endl;

	competition_file << "# Dispersal length:" << std::endl;
	for(i=0;i<num_species;i++) {
		competition_file << " " << dispersal_length[i];
	}
	competition_file << std::endl;

	competition_file << "# Intrinsic fecundity:" << std::endl;
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

	competition_file << "# fecundity transitivity row sum:" << std::endl;
	for(i=0;i<num_species;i++) {
		competition_file << " " << fecundity_row_sum[i];
	}
	competition_file << std::endl;

	competition_file << "# growth row sum:" << std::endl;
	for(i=0;i<num_species;i++) {
		competition_file << " " << growth_row_sum[i];
	}
	competition_file << std::endl;

	competition_file << "# Fecundity transitivity:" << std::endl;
	for(i=0;i<num_species;i++) {
		for(j=0;j<num_species;j++) {

			if(fecundity_transitivity[i][j] < 0)
				competition_file << fecundity_transitivity[i][j];
			else
				competition_file << " " << fecundity_transitivity[i][j];

			if(j!=num_species-1)
				competition_file << ", ";

			}
			competition_file << std::endl;
		}

	competition_file << "# Growth transitivity:" << std::endl;
	for(i=0;i<num_species;i++) {
		for(j=0;j<num_species;j++) {

			if(growth_transitivity[i][j] < 0)
				competition_file << growth_transitivity[i][j];
			else
				competition_file << " " << growth_transitivity[i][j];

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

	if(!properties_file.is_open()) {
			if(id==0)
				fprintf(stderr, "Error, could not open properties file for for saving\n");
			MPI_Finalize();
			exit(0);
	}

	properties_file <<  "# Fecundity imbalance mean:" << std::endl  << fecundity_imbalance_mean <<  std::endl;
	properties_file <<  "# Growth imbalance mean:" << std::endl  << growth_imbalance_mean <<  std::endl;
	properties_file <<  "# Fecundity relative intransitivity:" << std::endl  << fecundity_relative_intransitivity <<  std::endl;
	properties_file <<  "# Growth relative intransitivity:" << std::endl  << growth_relative_intransitivity <<  std::endl;
	properties_file <<  "# Fecundity-growth cross correlation:" << std::endl  << fecundity_growth_correlation <<  std::endl;

	properties_file.close();

	return;

}


int Simulation::getBoxSize() {
	return box_size;
}


int Simulation::getSpecies() {
	return num_species;
}


int Simulation::getMaxTimeStep() {
	return max_time_step;
}


int Simulation::getSite(int i, int j) {
	return box[i][j];
}


int Simulation::getNextSite(int i, int j) {
	return next_box[i][j];
}


int Simulation::getRandom() {
	random_count++;
	return global_random_generator();

}

unsigned long long Simulation::getRandomCount() {
	return random_count;
}


double Simulation::getDispersal(int i, int j, int s) {
	return dispersal_box[i][j][s];
}


double Simulation::getNextDispersal(int i, int j, int s) {
	return next_dispersal_box[i][j][s];
}


void Simulation::setSite(int i, int j, int s) {
	box[i][j] = s;
	return;
}


void Simulation::addSite(int i, int j, int s) {
	box[i][j] += s;
	return;
}

void Simulation::setDispersal(int i, int j, int s, double t) {
	dispersal_box[i][j][s] = t;
	return;
}

void Simulation::addDispersal(int i, int j, int s, double t) {
	dispersal_box[i][j][s] += t;
	return;
}

void Simulation::discardRandom(unsigned long long n) {
	global_random_generator.discard(n);
	random_count+=n;
}

