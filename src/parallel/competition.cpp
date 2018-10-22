
	 /**********************************************************
	 *
	 *			D.S Jamieson and N.L Kinlock, 2018			
	 *
	 *	These methods for the Simulation class set up random 
	 *	competition matrices, and manipulate and calculated their
	 *	properties.
	 *
	 ***********************************************************/

#include "simulation.h"

void Simulation::initializeUniformCompetition(){

	int i, j;

	std::uniform_real_distribution<double> dist(competition_lower_bound, competition_upper_bound);

	for( i = 0 ; i < num_species ; i++ ){
		for( j = 0 ; j < num_species ; j++ ){

			competition_fecundity[i][j] = dist(generateRandom());
			competition_growth[i][j] = dist(generateRandom());

		}
	}

	return;

}


void Simulation::initializeTNormalCompetition() {

	int i, j;

	std::normal_distribution<double> dist(competition_mean, competition_sdev);

	for(i=0; i<num_species; i++) {
		for(j=0; j<num_species;j++) {
			
			// Rejection sampling to obey bounds
			while(1) {
				competition_fecundity[i][j] = dist(generateRandom());
				if(competition_lower_bound <= competition_fecundity[i][j] && competition_upper_bound >= competition_fecundity[i][j])
					break;
			}

			while(1) {
				competition_growth[i][j] = dist(generateRandom());
				if(competition_lower_bound <= competition_growth[i][j] && competition_upper_bound >= competition_growth[i][j])
					break;			
			}
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

	// Rejection sampling to get correlation near competition_correlation
	while(1) {

		for( i = 0 ; i < num_species ; i++ ){
			for( j = 0 ; j < num_species ; j++ ){

				b = (double) bdist(generateRandom());
				x = udist(generateRandom());
				y = udist(generateRandom());
				z= udist(generateRandom());

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

	// Rejection sampling to get correlation near competition_correlation
	while(1) {

		for( i = 0 ; i < num_species ; i++ ){
			for( j = 0 ; j < num_species ; j++ ){

				// Rejection sampling to obey bounds
				while(1) {

					b = (double) bdist(generateRandom());
					x = ndist(generateRandom());
					y = ndist(generateRandom());
					z = ndist(generateRandom());

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
	double fecundity_fecundity_diff_upper, fecundity_fecundity_diff_lower, growth_fecundity_diff_upper, growth_fecundity_diff_lower;
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

			fecundity_fecundity_diff_upper = competition_fecundity[i][j] - fecundity_mean ;
 			fecundity_fecundity_diff_lower = competition_fecundity[j][i] - fecundity_mean;
			growth_fecundity_diff_upper = competition_growth[i][j] - growth_mean;
			growth_fecundity_diff_lower = competition_growth[j][i] - growth_mean;

			if(  fecundity_fecundity_diff_upper*fecundity_fecundity_diff_lower > 0.   ) {
				if(bdist(generateRandom())) {
					competition_fecundity[i][j] = fecundity_mean - fecundity_fecundity_diff_upper;
				}
			}
			else {
				if(!bdist(generateRandom())) {
					competition_fecundity[i][j] = fecundity_mean - fecundity_fecundity_diff_upper;
				}
			}

			if(  growth_fecundity_diff_upper*growth_fecundity_diff_lower > 0.   ) {
				if(bdist(generateRandom())) {
					competition_growth[i][j] = growth_mean - growth_fecundity_diff_upper;
				}
			}
			else {
				if(!bdist(generateRandom())) {
					competition_growth[i][j] = growth_mean - growth_fecundity_diff_upper;
				}
			}
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
	double top_mean = 0.;
	double bottom_mean = 0.;
	double top_variance = 0.;
	double bottom_variance = 0.;
	double fecundity_row_mean = 0.;	
	double growth_row_mean = 0.;
	double fecundity_variance = 0.;
	double growth_variance = 0.;

	int *top_row_sum = new int[num_species];
	int *bottom_row_sum = new int[num_species];
	if(!top_row_sum || !bottom_row_sum ) {
		fprintf(stderr, "Error, unable to allocate memory for transitivity calculation\n");
		MPI_Finalize();
		exit(-1);
	}

	for( i = 0 ; i < num_species ; i++ ){
		fecundity_row_sum[i] = 0.;
		growth_row_sum[i] = 0.;
	}

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

	delete[] top_row_sum;
	delete[] bottom_row_sum;

	for(i=0; i<num_species; i++) {
		for(j=i+1; j<num_species; j++) {

			if( competition_fecundity[i][j] < competition_fecundity[j][i]) {
				fecundity_transitivity[i][j] = 1.;
				fecundity_transitivity[j][i] = 0.;
			}
			else if( competition_fecundity[i][j] > competition_fecundity[j][i]  ) {
				fecundity_transitivity[i][j] = 0.;
				fecundity_transitivity[j][i] = 1.;
			}
			fecundity_row_sum[i] += fecundity_transitivity[i][j];
			fecundity_row_sum[j] += fecundity_transitivity[j][i];

			if( competition_growth[i][j] < competition_growth[j][i]) {
				growth_transitivity[i][j] = 1.;
				growth_transitivity[j][i] = 0.;
			}
			else if( competition_growth[i][j] > competition_growth[j][i]  ) {
				growth_transitivity[i][j] = 0.;
				growth_transitivity[j][i] = 1.;
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

	int i, j, k;
	double pivot = 0;

	for(i=0; i<num_species; i++) {
		 fecundity_row_sum[i] = 0.;
	}

	for(i=0; i<num_species; i++) {
		for(j=i+1; j<num_species; j++) {

			if( competition_fecundity[i][j] < competition_fecundity[j][i]) {
				fecundity_transitivity[i][j] = 1. ;
				fecundity_transitivity[j][i] = 0. ;
			}
			else if( competition_fecundity[i][j] > competition_fecundity[j][i]  ) {
				fecundity_transitivity[i][j] = 0.	;
				fecundity_transitivity[j][i] = 1.	;
			}

			fecundity_row_sum[i] += fecundity_transitivity[i][j];
			fecundity_row_sum[j] += fecundity_transitivity[j][i];

		
		}	
	}

	return;

}


void Simulation::getDiscreteGrowthTransitivity() {

	int i, j, k;
	double pivot = 0;

	for(i=0; i<num_species; i++) {
		 growth_row_sum[i] = 0.;
	}

	for(i=0; i<num_species; i++) {
		for(j=i+1; j<num_species; j++) {

			if( competition_growth[i][j] < competition_growth[j][i]) {
				growth_transitivity[i][j] = 1. ;
				growth_transitivity[j][i] = 0. ;
			}
			else if( competition_growth[i][j] > competition_growth[j][i]  ) {
				growth_transitivity[i][j] = 0.	;
				growth_transitivity[j][i] = 1.	;
			}
			growth_row_sum[i] += growth_transitivity[i][j];
			growth_row_sum[j] += growth_transitivity[j][i];

		}	
	}

	return;

}


void Simulation::setCompetitionTransitivity() {

	int i, j, k;
	int ih, jh;
	double hold;
	int *desired_fecundity_row_sum = new int[num_species];
	if( !desired_fecundity_row_sum ) {
		fprintf(stderr, "Error, unable to allocate memory for initializing transitive competition\n");
		MPI_Finalize();
		exit(-1);
	}
	if(fecundity_transitivity_type == 1.) {
		for(i = 0; i<num_species; i++)
			desired_fecundity_row_sum[i] = num_species - (i+1) ;
	}
	else if(fecundity_transitivity_type == -1.)  {
		if( num_species%2==0 ) {
			for(i = 0; i<num_species/2; i++) {
				desired_fecundity_row_sum[i] = num_species/2 ;
				desired_fecundity_row_sum[num_species-(i+1)] = num_species/2 -1 ; 
			}
		}
		else {
			for(i = 0; i<num_species; i++) {
				desired_fecundity_row_sum[i] = (num_species-1)/2 ;
			}
		}
	}
	else if( fecundity_transitivity_type == 0 ) {

		int * fecundity_hierarchy = new int[num_species];
		if(!fecundity_hierarchy ) {
			fprintf(stderr, "error, unable to allocate memory for setting transitivity hierarchy\n");
			MPI_Finalize();
			exit(-1);
		}
		for(i=0; i<num_species; i++) {
			fecundity_hierarchy[i] = 0;
		}

		setGrowthCompetitionTransitivity( fecundity_hierarchy );

		delete[] fecundity_hierarchy;

		return;

	}

	getDiscreteFecundityTransitivity();

	int * fecundity_hierarchy = new int[num_species];
	int * fecundity_diff = new int[num_species];
	if(!fecundity_hierarchy || !fecundity_diff ) {
		fprintf(stderr, "error, unable to allocate memory for setting transitivity hierarchy\n");
		MPI_Finalize();
		exit(-1);
	}

	for(i=0; i<num_species; i++)
		fecundity_hierarchy[i] = i;
	// Shuffle and calculate fecundity_difference from desired row sum
	for(i=num_species-1; i>=0 ;i--) {

		j = getRandom() % (i+1);
		k = fecundity_hierarchy[i];
		fecundity_hierarchy[i] = fecundity_hierarchy[j];
		fecundity_hierarchy[j] = k;
		ih = fecundity_hierarchy[i];
		fecundity_diff[ih] = (int) fecundity_row_sum[ih] - desired_fecundity_row_sum[i];
	
	} 

	// Fix row sums to desired value
	for(i=0; i<num_species; i++) {
	
		ih = fecundity_hierarchy[i];
		
		for(j=num_species-1; j>i; j--) {

			jh = fecundity_hierarchy[j];

			if( fecundity_diff[ih] < 0  &&  fecundity_transitivity[ih][jh] == 0   ) {
				hold = competition_fecundity[jh][ih];
				competition_fecundity[jh][ih] = competition_fecundity[ih][jh];
				competition_fecundity[ih][jh] = hold;
				fecundity_diff[ih]++;
				fecundity_diff[jh]--;
			}
			else if(fecundity_diff[ih] > 0  &&  fecundity_transitivity[ih][jh] == 1 ) {
				hold = competition_fecundity[jh][ih];
				competition_fecundity[jh][ih] = competition_fecundity[ih][jh];
				competition_fecundity[ih][jh] = hold;
				fecundity_diff[ih]--;
				fecundity_diff[jh]++;
			}

		}

	}

	getDiscreteFecundityTransitivity();

	// In intransitive case, one transitive loop could remain, this checks it
	ih = -1;
	jh = -1;	
	for(i=0; i<num_species; i++) {
		ih = fecundity_hierarchy[i];
		fecundity_diff[ih] = (int) fecundity_row_sum[ih] - desired_fecundity_row_sum[i];
		if(fecundity_diff[ih] == 1)
			break;
	}  

	for(i=0; i<num_species; i++) {
		jh = fecundity_hierarchy[i];
		fecundity_diff[jh] = (int) fecundity_row_sum[jh] - desired_fecundity_row_sum[i];
		if(fecundity_diff[jh] == -1)
			break;
	}

	// If transitive triplet remains, this fixes it
	if( ih != -1  ) {
		for(i==0; i<num_species; i++) {

			if( i == ih || i == jh )
				continue;

			if( fecundity_transitivity[ih][i] == 1 &&  fecundity_transitivity[jh][i] == 0 ) {
				hold = competition_fecundity[ih][i];
				competition_fecundity[ih][i] = competition_fecundity[i][ih];
				competition_fecundity[i][ih] = hold;
				hold = competition_fecundity[jh][i];
				competition_fecundity[jh][i] = competition_fecundity[i][jh];
				competition_fecundity[i][jh] = hold;
				break;

			}
		}
	}


	delete[] desired_fecundity_row_sum;
	delete[] fecundity_diff;

	getDiscreteTransitivity();

	// Fix growth competition
	if( growth_transitivity_type == fecundity_transitivity_type && fecundity_growth_relative_hierarchy != 0.  ) {

		for(i=0; i<num_species; i++) {
			for(j=i+1; j<num_species; j++) {

				if( i == j )
					continue;

				if( fecundity_transitivity[i][j] != growth_transitivity[i][j] && fecundity_growth_relative_hierarchy == 1  ) {
					hold = competition_growth[j][i];
					competition_growth[j][i] = competition_growth[i][j];
					competition_growth[i][j] = hold;
				}
				else if( fecundity_transitivity[i][j] == growth_transitivity[i][j] && fecundity_growth_relative_hierarchy == -1  ) {
					hold = competition_growth[j][i];
					competition_growth[j][i] = competition_growth[i][j];
					competition_growth[i][j] = hold;
				}

			}
		}
	}
	else {

		setGrowthCompetitionTransitivity( fecundity_hierarchy );

	}

	delete[] fecundity_hierarchy;


	return;

}


void Simulation::setGrowthCompetitionTransitivity( int * fecundity_hierarchy ) {

	int i, j, k;
	int ih, jh;
	double hold;
	int *desired_growth_row_sum = new int[num_species];

	if( !desired_growth_row_sum ) {
		fprintf(stderr, "Error, unable to allocate memory for initializing transitive competition\n");
		MPI_Finalize();
		exit(-1);
	}
	if(growth_transitivity_type == 1.){
		for(i = 0; i<num_species; i++)
			desired_growth_row_sum[i] = num_species - (i+1) ;
	}
	else if( growth_transitivity_type == -1. )  {
		if( num_species%2==0 ) {
			for(i = 0; i<num_species/2; i++) {
				desired_growth_row_sum[i] = num_species/2 ;
				desired_growth_row_sum[num_species-(i+1)] = num_species/2 -1 ; 
			}
		}
		else {
			for(i = 0; i<num_species; i++) {
				desired_growth_row_sum[i] = (num_species-1)/2 ;
			}
		}
	}

	getDiscreteGrowthTransitivity();

	int * growth_hierarchy = new int[num_species];
	int * growth_diff = new int[num_species];
	if(!growth_hierarchy || !growth_diff ) {
		fprintf(stderr, "Error, unable to allocate memory for setting transitivity hierarchy\n");
		MPI_Finalize();
		exit(-1);
	}

	if( fecundity_growth_relative_hierarchy == 0) {

		for(i=0; i<num_species; i++)
			growth_hierarchy[i] = i;
		// Shuffle and calculate growth_difference from desired row sum
		for(i=num_species-1; i>=0 ;i--) {

			j = getRandom() % (i+1);
			k = growth_hierarchy[i];
			growth_hierarchy[i] = growth_hierarchy[j];
			growth_hierarchy[j] = k;
			ih = growth_hierarchy[i];
			growth_diff[ih] = (int) growth_row_sum[ih] - desired_growth_row_sum[i];
		
		} 
	}
	else {

		for(i=0; i<num_species; i++) {

			if(  fecundity_growth_relative_hierarchy == 1  )
				growth_hierarchy[i] = fecundity_hierarchy[i];
			if(  fecundity_growth_relative_hierarchy == -1  )
				growth_hierarchy[i] = fecundity_hierarchy[num_species - (i+1)];

			ih = growth_hierarchy[i];
			growth_diff[ih] = (int) growth_row_sum[ih] - desired_growth_row_sum[i];

		}
	}

	// Fix row sums to desired value
	for(i=0; i<num_species; i++) {
	
		ih = growth_hierarchy[i];
		
		for(j=num_species-1; j>i; j--) {

			jh = growth_hierarchy[j];

			if( growth_diff[ih] < 0  &&  growth_transitivity[ih][jh] == 0   ) {
				hold = competition_growth[jh][ih];
				competition_growth[jh][ih] = competition_growth[ih][jh];
				competition_growth[ih][jh] = hold;
				growth_diff[ih]++;
				growth_diff[jh]--;
			}
			else if(growth_diff[ih] > 0  &&  growth_transitivity[ih][jh] == 1 ) {
				hold = competition_growth[jh][ih];
				competition_growth[jh][ih] = competition_growth[ih][jh];
				competition_growth[ih][jh] = hold;
				growth_diff[ih]--;
				growth_diff[jh]++;
			}

		}

	}

	getDiscreteGrowthTransitivity();

	// In intransitive case, one transitive loop could remain, this checks it
	ih = -1;
	jh = -1;	
	for(i=0; i<num_species; i++) {
		ih = growth_hierarchy[i];
		growth_diff[ih] = (int) growth_row_sum[ih] - desired_growth_row_sum[i];
		if(growth_diff[ih] == 1)
			break;
	}  

	for(i=0; i<num_species; i++) {
		jh = growth_hierarchy[i];
		growth_diff[jh] = (int) growth_row_sum[jh] - desired_growth_row_sum[i];
		if(growth_diff[jh] == -1)
			break;
	}

	// If transitive triplet remains, this fixes it
	if( ih != -1  ) {
		for(i==0; i<num_species; i++) {

			if( i == ih || i == jh )
				continue;

			if( growth_transitivity[ih][i] == 1 &&  growth_transitivity[jh][i] == 0 ) {
				hold = competition_growth[ih][i];
				competition_growth[ih][i] = competition_growth[i][ih];
				competition_growth[i][ih] = hold;
				hold = competition_growth[jh][i];
				competition_growth[jh][i] = competition_growth[i][jh];
				competition_growth[i][jh] = hold;
				break;

			}
		}
	}


	delete[] desired_growth_row_sum;
	delete[] growth_hierarchy;
	delete[] growth_diff;

	return;

}

