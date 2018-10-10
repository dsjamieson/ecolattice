// Percentage juvenile in initial conditions option
// optaional condiiton function for parameter file
// initialization defaults
// work out overfilling box problem

#include <cmath>
#include <random>
#include "simulation.h"

void Simulation::initializeUniformCompetition(){

	int i, j;

	// Initialize distributions
	std::uniform_real_distribution<double> dist(competition_lower_bound, competition_upper_bound);

	// Allocate arrays
	competition_fecundity = new double *[num_species];
	for( i = 0 ; i < num_species ; i++ ){
		competition_fecundity[i] = new double[num_species];
		for(j=0;j<num_species;j++) competition_fecundity[i][j] = 0.;
	}
	competition_growth = new double *[num_species];
	for( i = 0 ; i < num_species ; i++ ){
		competition_growth[i] = new double[num_species];
		for(j=0;j<num_species;j++) competition_growth[i][j] = 0.;
	}

	for( i = 0 ; i < num_species ; i++ ){
		for( j = 0 ; j < num_species ; j++ ){

			competition_fecundity[i][j] = dist(global_random_generator);
			competition_growth[i][j] = dist(global_random_generator);

		}
	}

	return;

}

void Simulation::initializeTNormalCompetition() {

	int i, j;
	double mean, comphold;

	// Allocate arrays
	competition_fecundity = new double *[num_species];
	for( i = 0 ; i < num_species ; i++ ){
		competition_fecundity[i] = new double[num_species];
		for(j=0;j<num_species;j++) competition_fecundity[i][j] = 0.;
	}
	competition_growth = new double *[num_species];
	for( i = 0 ; i < num_species ; i++ ){
		competition_growth[i] = new double[num_species];
		for(j=0;j<num_species;j++) competition_growth[i][j] = 0.;
	}

	mean = (competition_upper_bound + competition_lower_bound)/2.;
	std::normal_distribution<double> dist(mean, competition_sdev);


	for(i=0; i<num_species; i++) {
		for(j=0; j<num_species;j++) {

			while(1) {
				comphold = dist(global_random_generator);
				if(competition_lower_bound <= comphold && competition_upper_bound >= comphold) {
						break;
					}
				}
			competition_fecundity[i][j] = comphold;

			while(1) {
				comphold = dist(global_random_generator);
				if(competition_lower_bound <= comphold && competition_upper_bound >= comphold) {
						break;
					}
				}
			competition_growth[i][j] = comphold;

		}
	}

	return;

}




void Simulation::initializeUniformCorrelatedCompetition(){

	int i, j;
	double b, x, y, z, umean;

	umean = (competition_lower_bound + competition_upper_bound)/2.;
	std::uniform_real_distribution<double> udist(competition_lower_bound, competition_upper_bound);
	std::bernoulli_distribution bdist(fabs(competition_correlation));

	// Allocate arrays
	competition_fecundity = new double *[num_species];
	for( i = 0 ; i < num_species ; i++ ){
		competition_fecundity[i] = new double[num_species];
		for(j=0;j<num_species;j++) competition_fecundity[i][j] = 0.;
	}
	competition_growth = new double *[num_species];
	for( i = 0 ; i < num_species ; i++ ){
		competition_growth[i] = new double[num_species];
		for(j=0;j<num_species;j++) competition_growth[i][j] = 0.;
	}

	while(1) {


		for( i = 0 ; i < num_species ; i++ ){
			for( j = 0 ; j < num_species ; j++ ){

				b = (double) bdist(global_random_generator);
				x = udist(global_random_generator);
				y = udist(global_random_generator);
				z= udist(global_random_generator);

				if(competition_correlation > 0 ) {
					competition_fecundity[i][j] = b*x + (1. - b)*y;
					competition_growth[i][j] = b*x + (1.-b)*z;
				}
				else {
					competition_fecundity[i][j] = b*x + (1. - b)*y;
					competition_growth[i][j] = b*(-x+2.*umean) + (1.-b)*z;
				}
			}
		}

	getFecundityGrowthCorrelation();
	if(fabs(fecundity_growth_correlation-competition_correlation) <= 0.05)
		break;
	}
	
	return;
}

void Simulation::initializeTNormalCorrelatedCompetition(){

	int i, j;
	double b, x, y, z, nmean;

	nmean = (competition_lower_bound + competition_upper_bound)/2.;
	std::normal_distribution<double> ndist(nmean, competition_sdev);
	std::bernoulli_distribution bdist(fabs(competition_correlation));

	// Allocate arrays
	competition_fecundity = new double *[num_species];
	for( i = 0 ; i < num_species ; i++ ){
		competition_fecundity[i] = new double[num_species];
		for(j=0;j<num_species;j++) competition_fecundity[i][j] = 0.;
	}
	competition_growth = new double *[num_species];
	for( i = 0 ; i < num_species ; i++ ){
		competition_growth[i] = new double[num_species];
		for(j=0;j<num_species;j++) competition_growth[i][j] = 0.;
	}

	while(1) {

		for( i = 0 ; i < num_species ; i++ ){
			for( j = 0 ; j < num_species ; j++ ){

				while(1) {

					b = (double) bdist(global_random_generator);
					x = ndist(global_random_generator);
					y = ndist(global_random_generator);
					z = ndist(global_random_generator);

					if(competition_correlation > 0 ) {
						competition_fecundity[i][j] = b*x + (1. - b)*y;
						competition_growth[i][j] = b*x + (1.- b)*z;
					}
					else {
						competition_fecundity[i][j] = b*x + (1. - b)*y;
						competition_growth[i][j] = b*(-x + 2.*nmean) + (1.- b)*z;
					}
				

					if(competition_lower_bound <= competition_fecundity[i][j] && competition_upper_bound >= competition_fecundity[i][j]) {
						if(competition_lower_bound <= competition_growth[i][j] && competition_upper_bound >= competition_growth[i][j] ) {
							break;
						}
					}
				}



			}
		}

		getFecundityGrowthCorrelation();
		if(fabs(fecundity_growth_correlation-competition_correlation) <= 0.05)
			break;
	}

	return;

}

void Simulation::imbalanceCompetition() {

	int i, j;
	double fecundity_diff_upper, fecundity_diff_lower, growth_diff_upper, growth_diff_lower;
	double fecundity_mean = 0;
	double growth_mean = 0;
	double pre_check = 0.;
	double imbalance_check = 0.;
	std::bernoulli_distribution bdist(imbalance);

	for(i=0;i<num_species;i++) {
		for(j=i+1;j<num_species;j++) {
			if(i!=j){
				fecundity_mean+=competition_fecundity[i][j];
				growth_mean+=competition_growth[i][j];
			}
		}

	}
	fecundity_mean = fecundity_mean/num_species/(num_species-1.);
	growth_mean = growth_mean/num_species/(num_species-1.);

	for(i=0;i<num_species;i++) {
		for(j=i+1;j<num_species;j++) {

			fecundity_diff_upper = competition_fecundity[i][j] - fecundity_mean ;
 			fecundity_diff_lower = competition_fecundity[j][i] - fecundity_mean;
			growth_diff_upper = competition_growth[i][j] - growth_mean;
			growth_diff_lower = competition_growth[j][i] - growth_mean;

			if(  fecundity_diff_upper*fecundity_diff_lower > 0.   ) {
				if(bdist(global_random_generator)) {
					competition_fecundity[i][j] = fecundity_mean - fecundity_diff_upper;
				}
			}
			else {
				if(!bdist(global_random_generator)) {
					competition_fecundity[i][j] = fecundity_mean - fecundity_diff_upper;
				}
			}

			if(  growth_diff_upper*growth_diff_lower > 0.   ) {
				if(bdist(global_random_generator)) {
					competition_growth[i][j] = growth_mean - growth_diff_upper;
				}
			}
			else {
				if(!bdist(global_random_generator)) {
					competition_growth[i][j] = growth_mean - growth_diff_upper;
				}
			}
		}
	}

	return;

}

void Simulation::getImbalanceMean() {
	
	int i, j;
	fecundity_imbalance_mean = 0.;
	growth_imbalance_mean = 0.;

	for(i=0;i<num_species;i++) {
		for(j=i+1;j<num_species;j++) {
			if(i!=j){
				fecundity_imbalance_mean += fabs(competition_fecundity[i][j] - competition_fecundity[j][i]);
				growth_imbalance_mean += fabs(competition_growth[i][j] - competition_growth[j][i]);
			}
		}
	}
	fecundity_imbalance_mean = fecundity_imbalance_mean/num_species/(num_species-1.);
	growth_imbalance_mean = growth_imbalance_mean/num_species/(num_species-1.);

	return;

}

void Simulation::getDiscreteTransitivity() {

	int i, j;
	int top_row_sum[num_species];
	int	bottom_row_sum[num_species];
	double top_mean = 0.;
	double bottom_mean = 0.;
	double top_variance = 0.;
	double bottom_variance = 0.;
	double fecundity_row_mean = 0.;	
	double growth_row_mean = 0.;
	double fecundity_variance = 0.;
	double growth_variance = 0.;

	for(i=0;i<num_species;i++) {
		top_row_sum[i] = num_species - (i+1);
		if(num_species%2==0) {
			if(i<num_species/2)
				bottom_row_sum[i] = num_species/2;
			else
				bottom_row_sum[i] = num_species/2. - 1;
		}
		else
			bottom_row_sum[i] = (num_species-1)/2;

		top_mean += top_row_sum[i];
		bottom_mean += bottom_row_sum[i];
	}
	top_mean = top_mean/num_species;
	bottom_mean = bottom_mean/num_species;

	for(i=0;i<num_species;i++) {
		top_variance += pow(top_row_sum[i]-top_mean, 2);
		bottom_variance += pow(bottom_row_sum[i]-bottom_mean, 2);
	}
	top_variance = top_variance/(num_species-1);
	bottom_variance = bottom_variance/(num_species-1);


	fecundity_row_sum = new double[num_species];
	growth_row_sum = new double[num_species];
	fecundity_transitivity = new double *[num_species];
	growth_transitivity = new double *[num_species];
	for( i = 0 ; i < num_species ; i++ ){
		fecundity_row_sum[i] = 0.;
		growth_row_sum[i] = 0.;
		fecundity_transitivity[i] = new double[num_species];
		growth_transitivity[i] = new double[num_species];
		for(j=0;j<num_species;j++) {
			fecundity_transitivity[i][j] = 0.;
			growth_transitivity[i][j] = 0.;
		}
	}

	for(i=0; i<num_species; i++) {
		for(j=i+1; j<num_species; j++) {

			if( competition_fecundity[i][j] < competition_fecundity[j][i]) {
				fecundity_transitivity[i][j] = 1. ;
				fecundity_transitivity[j][i] = 0 ;
			}
			else if( competition_fecundity[i][j] > competition_fecundity[j][i]  ) {
				fecundity_transitivity[i][j] = 0	;
				fecundity_transitivity[j][i] = 1.	;
			}
			fecundity_row_sum[i] += fecundity_transitivity[i][j];
			fecundity_row_sum[j] += fecundity_transitivity[j][i];

			if( competition_growth[i][j] < competition_growth[j][i]) {
				growth_transitivity[i][j] = 1. ;
				growth_transitivity[j][i] = 0. ;
			}
			else if( competition_growth[i][j] > competition_growth[j][i]  ) {
				growth_transitivity[i][j] = 0	;
				growth_transitivity[j][i] = 1.	;
			}
			growth_row_sum[i] += growth_transitivity[i][j];
			growth_row_sum[j] += growth_transitivity[j][i];

		}	
	}

	for(i=0;i<num_species;i++) {
		fecundity_row_mean += fecundity_row_sum[i];
		growth_row_mean += growth_row_sum[i];
	}
	fecundity_row_mean = fecundity_row_mean/num_species;
	growth_row_mean = growth_row_mean/num_species;

	for(i=0;i<num_species;i++) {
		fecundity_variance += pow(fecundity_row_sum[i]-fecundity_row_mean, 2);
		growth_variance += pow(growth_row_sum[i]-growth_row_mean, 2);
	}
	fecundity_variance = fecundity_variance/(num_species-1);
	growth_variance = growth_variance/(num_species-1);

	fecundity_relative_intransitivity = 1. - (fecundity_variance-bottom_variance)/(top_variance-bottom_variance);
	growth_relative_intransitivity = 1. -  (growth_variance-bottom_variance)/(top_variance-bottom_variance);

	return;

}

void Simulation::getContinuousTransitivity() {

	int i, j;

	fecundity_row_sum = new double[num_species];
	growth_row_sum = new double[num_species];
	fecundity_transitivity = new double *[num_species];
	growth_transitivity = new double *[num_species];
	for( i = 0 ; i < num_species ; i++ ){
		fecundity_row_sum[i] = 0.;
		growth_row_sum[i] = 0.;
		fecundity_transitivity[i] = new double[num_species];
		growth_transitivity[i] = new double[num_species];
		for(j=0;j<num_species;j++) {
			fecundity_transitivity[i][j] = 0.;
			growth_transitivity[i][j] = 0.;
		}
	}

	for(i=0; i<num_species; i++) {
		for(j=i+1; j<num_species; j++) {

			fecundity_transitivity[i][j] = (- competition_fecundity[i][j] + competition_fecundity[j][i] )
											/(competition_upper_bound-competition_lower_bound);
			fecundity_transitivity[j][i] = -fecundity_transitivity[i][j];
			fecundity_row_sum[i] += fecundity_transitivity[i][j];
			fecundity_row_sum[j] += fecundity_transitivity[j][i];

			growth_transitivity[i][j] = ( - competition_growth[i][j] + competition_growth[j][i] )
											/ ( competition_upper_bound-competition_lower_bound  );

			growth_transitivity[j][i] = -growth_transitivity[i][j];
			growth_row_sum[i] += growth_transitivity[i][j];
			growth_row_sum[j] += growth_transitivity[j][i];

		}	
	}

	return;

}

void Simulation::getDiscreteFecundityTransitivity() {



	int i, j;

	fecundity_row_sum = new double[num_species];
	fecundity_transitivity = new double *[num_species];
	for( i = 0 ; i < num_species ; i++ ){
		fecundity_row_sum[i] = 0.;
		fecundity_transitivity[i] = new double[num_species];
		for(j=0;j<num_species;j++) {
			fecundity_transitivity[i][j] = 0.;
		}
	}

	for(i=0; i<num_species; i++) {
		for(j=i+1; j<num_species; j++) {

			if( competition_fecundity[i][j] < competition_fecundity[j][i]) {
				fecundity_transitivity[i][j] = 1. ;
				fecundity_transitivity[j][i] = 0 ;
			}
			else if( competition_fecundity[i][j] > competition_fecundity[j][i]  ) {
				fecundity_transitivity[i][j] = 0	;
				fecundity_transitivity[j][i] = 1.	;
			}
			fecundity_row_sum[i] += fecundity_transitivity[i][j];
			fecundity_row_sum[j] += fecundity_transitivity[j][i];
		}	
	}

	double pivot = 0;
	for(i=0; i<num_species; i++) {
		for(j=i+1; j<num_species; j++) {
			if(fecundity_row_sum[j] > fecundity_row_sum[i]){
				pivot = fecundity_row_sum[i];
				fecundity_row_sum[i] = fecundity_row_sum[j];
				fecundity_row_sum[j] = pivot; 
			}
		}
	}


	return;

}

void Simulation::initializeUniformTransitiveCompetition() {

	int i, j;
	int set_row_sum[num_species];
	int top_row_sum[num_species];
	int	bottom_row_sum[num_species];

	// Allocate arrays
	competition_fecundity = new double *[num_species];
	for( i = 0 ; i < num_species ; i++ ){
		competition_fecundity[i] = new double[num_species];
		for(j=0;j<num_species;j++) competition_fecundity[i][j] = 0.;
	}
	competition_growth = new double *[num_species];
	for( i = 0 ; i < num_species ; i++ ){
		competition_growth[i] = new double[num_species];
		for(j=0;j<num_species;j++) competition_growth[i][j] = 0.;
	}

	if(transitivity == 1){
		for(i = 0; i<num_species; i++)
			set_row_sum[i] = num_species - (i+1) ;
	}
	else{
		if( num_species%2==0 ) {
			for(i = 0; i<num_species/2; i++) {
				set_row_sum[i] = num_species/2 ;
				set_row_sum[num_species-(i+1)] = num_species/2 -1 ; 
			}
		}
		else {
			for(i = 0; i<num_species; i++) {
				set_row_sum[i] = (num_species-1)/2 ;
			}
		}
	}

	std::uniform_real_distribution<double> dist(competition_lower_bound, competition_upper_bound);

	while(1) {

		for(i=0; i<num_species; i++) {
			for(j=0; j<num_species;j++) {
				competition_fecundity[i][j] = dist(global_random_generator);
			}
		}

		getDiscreteFecundityTransitivity();

		j = 1;
		for(i=0; i<num_species; i++) {
			if( set_row_sum[i] != fecundity_row_sum[i]) {
				j=0;
				break;
			}
		}

		if(j==1)
			break;
	
	}

	for(i=0; i<num_species; i++) {
		for(j=0; j<num_species;j++) {
			competition_growth[i][j] = competition_fecundity[i][j];
		}
	}

	return;

}


void Simulation::initializeTNormalTransitiveCompetition() {

	int i, j;
	int set_row_sums[num_species];
	double mean, comphold;

	// Allocate arrays
	competition_fecundity = new double *[num_species];
	for( i = 0 ; i < num_species ; i++ ){
		competition_fecundity[i] = new double[num_species];
		for(j=0;j<num_species;j++) competition_fecundity[i][j] = 0.;
	}
	competition_growth = new double *[num_species];
	for( i = 0 ; i < num_species ; i++ ){
		competition_growth[i] = new double[num_species];
		for(j=0;j<num_species;j++) competition_growth[i][j] = 0.;
	}

	if(transitivity == 1){
		for(i = 0; i<num_species; i++)
			set_row_sums[i] = num_species - (i+1) ;
	}
	else{
		if( num_species%2==0 ) {
			for(i = 0; i<num_species/2; i++) {
				set_row_sums[i] = num_species/2 ;
				set_row_sums[num_species-(i+1)] = num_species/2 -1 ; 
			}
		}
		else {
			for(i = 0; i<num_species; i++) {
				set_row_sums[i] = (num_species-1)/2 ;
			}
		}
	}

	mean = (competition_upper_bound + competition_lower_bound)/2.;
	std::normal_distribution<double> dist(mean, competition_sdev);

	while(1) {

		for(i=0; i<num_species; i++) {
			for(j=0; j<num_species;j++) {
				while(1) {
					comphold = dist(global_random_generator);
					if(competition_lower_bound <= comphold && competition_upper_bound >= comphold) {
							break;
						}
					}
				competition_fecundity[i][j] = comphold;
			}
		}

		getDiscreteFecundityTransitivity();

		j = 1;
		for(i=0; i<num_species; i++) {
			if( set_row_sums[i] != fecundity_row_sum[i]) {
				j=0;
				break;
			}
		}

		if(j==1)
			break;
	
	}

	for(i=0; i<num_species; i++) {
		for(j=0; j<num_species;j++) {		
			competition_growth[i][j] = competition_fecundity[i][j];
		}
	}

	return;



}

void Simulation::getFecundityGrowthCorrelation() {

	int i, j;
	double covariance = 0.;
	double fecundity_variance = 0.;
	double growth_variance = 0.;
	double fecundity_mean = 0.;
	double growth_mean = 0.;

	for(i=0;i<num_species;i++) {
		for(j=0;j<num_species;j++) {
			fecundity_mean += competition_fecundity[i][j];
			growth_mean += competition_growth[i][j];
		}
	}
	fecundity_mean = fecundity_mean/num_species/num_species;
	growth_mean = growth_mean/num_species/num_species;

	for(i=0;i<num_species;i++) {
		for(j=0;j<num_species;j++) {
			covariance += (competition_fecundity[i][j]-fecundity_mean)*(competition_growth[i][j]-growth_mean);
			fecundity_variance += pow(competition_fecundity[i][j]-fecundity_mean, 2);
			growth_variance += pow(competition_growth[i][j]-growth_mean,2);
		}
	}
	covariance = covariance/(num_species*num_species-1.);
	fecundity_variance = fecundity_variance/(num_species*num_species-1.);
	growth_variance = growth_variance/(num_species*num_species-1.);

	fecundity_growth_correlation = covariance/sqrt(fecundity_variance*growth_variance);




	return;

}

