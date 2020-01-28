
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2018			
	 *
	 *		methods for the Simulation class. set up random 
	 *		growth and fecundity competition matrices; 
	 *		manipulate imbalance, transitivity, and correlation;
 	 *		and calculate properties.
	 *
	 ***********************************************************/

#include "simulation.h"

void Simulation::initializeUniformCompetition(){
	/* set up growth (competition_growth) and fecundity (competition_fecundity) 
	competition matrices as random draws from a uniform distribution */

	int i, j;

	// random draw from Uniform(competition_lower_bound, competition_upper_bound)
	// separate draws for growth and fecundity competition matrices
	for (i = 0; i < num_species ; i++){
		for (j = 0; j < num_species; j++){
			if (i == j) {
			competition_fecundity[i][j] = getRandomUniformReal(competition_diag_lower_bound, competition_diag_upper_bound);
			competition_growth[i][j] = getRandomUniformReal(competition_diag_lower_bound, competition_diag_upper_bound);
			}
			else {
			competition_fecundity[i][j] = getRandomUniformReal(competition_lower_bound, competition_upper_bound);
			competition_growth[i][j] = getRandomUniformReal(competition_lower_bound, competition_upper_bound);
			}
		}
	}
	return;
}


void Simulation::initializeTNormalCompetition() {
	/* set up growth (competition_growth) and fecundity (competition_fecundity) 
	competition matrices as random draws from a truncated normal distribution */

	int i, j;

	// random draw from truncated Normal(competition_mean, competition_sdev) with bounds [competition_lower_bound, competition_upper_bound]
	// separate draws for growth and fecundity competition matrices
	for (i = 0; i < num_species; i++) {
		for (j = 0; j < num_species; j++) {
			
			// rejection sampling to obey bounds
			while (1) {
				if (i == j) {
				competition_fecundity[i][j] = getRandomNormal(competition_diag_mean, competition_diag_sdev);
				}
				else {
				competition_fecundity[i][j] = getRandomNormal(competition_mean, competition_sdev);
				}
				if(competition_lower_bound <= competition_fecundity[i][j] && competition_upper_bound >= competition_fecundity[i][j])
					break;
			}

			while (1) {
				if (i == j) {
				competition_growth[i][j] = getRandomNormal(competition_diag_mean, competition_diag_sdev);
				}
				else {
				competition_growth[i][j] = getRandomNormal(competition_mean, competition_sdev);
				}
				if (competition_lower_bound <= competition_growth[i][j] && competition_upper_bound >= competition_growth[i][j])
					break;			
			}
		}
	}
	return;
}


void Simulation::initializeUniformCorrelatedCompetition(){
	/* set up growth (competition_growth) and fecundity (competition_fecundity)
	competition matrices as random draws from a uniform distribution,
	and manipulate matrices so that they are correlated or anticorrelated.
	this represents a trade-off between growth and reproduction. */

	int i, j;
	double b, x, y, z;

	std::bernoulli_distribution bdist(fabs(competition_correlation));

	// draw three random uniform vectors: x, y, z
	// correlation is the Bernoulli probability that the matching elements in the two matrices are identical
	while (1) {

		for (i = 0; i < num_species ; i++ ) {
			for (j = 0; j < num_species; j++ ) {
				b = (double) bdist(generateRandom());
				if (i == j) {
					x = getRandomUniformReal(competition_diag_lower_bound, competition_diag_upper_bound);
					y = getRandomUniformReal(competition_diag_lower_bound, competition_diag_upper_bound);
					z = getRandomUniformReal(competition_diag_lower_bound, competition_diag_upper_bound);
				}
				else {
					x = getRandomUniformReal(competition_lower_bound, competition_upper_bound);
					y = getRandomUniformReal(competition_lower_bound, competition_upper_bound);
					z = getRandomUniformReal(competition_lower_bound, competition_upper_bound);
				}
				if (competition_correlation > 0 ) {
					competition_fecundity[i][j] = b * x + (1. - b) * y;
					competition_growth[i][j] = b * x + (1. - b) * z;
				}
				else {
					competition_fecundity[i][j] = b * x + (1. - b) * y;
					competition_growth[i][j] = b * (-x + 2. * competition_mean) + (1. - b) * z;
				}
			}
		}
	getFecundityGrowthCorrelation(); // calculate correlation of the two matrices

	// rejection sampling to get correlation near competition_correlation
	// control threshold for accepted correlation, currently 0.05
	if (fabs(fecundity_growth_correlation - competition_correlation) <= 0.05)
		break;
	}	
	return;
}


void Simulation::initializeTNormalCorrelatedCompetition(){
	/* set up growth (competition_growth) and fecundity (competition_fecundity)
	competition matrices as random draws from a truncated normal
	distribution, and manipulate matrices so that they are correlated or 
	anticorrelated. this represents a trade-off between growth and reproduction. */

	int i, j;
	double b, x, y, z;

	std::bernoulli_distribution bdist(fabs(competition_correlation));

	// draw three random uniform vectors: x, y, z
	// correlation is the Bernoulli probability that the matching elements in the two matrices are identical
	while(1) {

		for (i = 0; i < num_species; i++) {
			for (j = 0; j < num_species; j++) {

				while (1) {
					b = (double) bdist(generateRandom());
					if (i == j) {
						x = getRandomNormal(competition_diag_mean, competition_diag_sdev);
						y = getRandomNormal(competition_diag_mean, competition_diag_sdev);
						z = getRandomNormal(competition_diag_mean, competition_diag_sdev);
					}
					else {
						x = getRandomNormal(competition_mean, competition_sdev);
						y = getRandomNormal(competition_mean, competition_sdev);
						z = getRandomNormal(competition_mean, competition_sdev);
					}
					if( competition_correlation > 0 ) {
						competition_fecundity[i][j] = b * x + (1. - b) * y;
						competition_growth[i][j] = b * x + (1.- b) * z;
					}
					else {
						competition_fecundity[i][j] = b * x + (1. - b) * y;
						competition_growth[i][j] = b * (-x + 2. * competition_mean) + (1. - b) * z;
					}
				
					// first-level of rejection sampling to obey bounds of truncated normal
					if (i == j) {
						if (competition_diag_lower_bound <= competition_fecundity[i][j] && competition_diag_upper_bound >= competition_fecundity[i][j]) {
							if( competition_diag_lower_bound <= competition_growth[i][j] && competition_diag_upper_bound >= competition_growth[i][j] ) {
								break;
							}
						}
					}
					else {
						if (competition_lower_bound <= competition_fecundity[i][j] && competition_upper_bound >= competition_fecundity[i][j]) {
							if( competition_lower_bound <= competition_growth[i][j] && competition_upper_bound >= competition_growth[i][j] ) {
								break;
							}	
						}
					}
				}
			}
		}
		getFecundityGrowthCorrelation(); // calculate correlation of the two matrices

		// second-level of rejection sampling to get correlation near competition_correlation
		// control threshold for accepted correlation, currently 0.05
		if(fabs(fecundity_growth_correlation - competition_correlation) <= 0.05)
			break;
	}
	return;
}

void Simulation::imbalanceCompetition() {
	/* manipulate imbalance of (already initialized) growth and fecundity competition matrices.
	imbalance is governed by 'imbalance' parameter, the Bernoulli probability of increasing
	the imbalance of paired elements. imbalance refers to one species having a greater effect 
	on a second species than the second has on the first. */

	int i, j;
	double fecundity_fecundity_diff_upper, fecundity_fecundity_diff_lower, growth_fecundity_diff_upper, growth_fecundity_diff_lower;
	double fecundity_mean = 0;
	double growth_mean = 0;
	double pre_check = 0.;
	double imbalance_check = 0.;
	std::bernoulli_distribution bdist(imbalance);

	// calculate mean of growth and fecundity competition matrices
	for (i = 0; i < num_species; i++) {
		for (j = i + 1; j < num_species; j++) {
			if (i != j) {
				fecundity_mean += competition_fecundity[i][j];
				growth_mean += competition_growth[i][j];
			}
		}
	}
	fecundity_mean = fecundity_mean / num_species / (num_species - 1.);
	growth_mean = growth_mean / num_species / (num_species - 1.);

	for (i = 0; i < num_species; i++) {
		for (j = i + 1; j < num_species; j++) {
			// calculate difference of each element from mean of matrix
			fecundity_fecundity_diff_upper = competition_fecundity[i][j] - fecundity_mean;
 			fecundity_fecundity_diff_lower = competition_fecundity[j][i] - fecundity_mean;
			growth_fecundity_diff_upper = competition_growth[i][j] - growth_mean;
			growth_fecundity_diff_lower = competition_growth[j][i] - growth_mean;
			
			// if paired elements (i, j and j, i) are balanced (on same side of mean)
			// imbalance them with Bernoulli probability equal to imbalance parameter
			// how to imbalance: reflect upper triangular element to the opposite side of the mean
			if (fecundity_fecundity_diff_upper * fecundity_fecundity_diff_lower > 0.) {
				if (bdist(generateRandom())) {
					competition_fecundity[i][j] = fecundity_mean - fecundity_fecundity_diff_upper;
				}
			}

			// if paired elements are already imbalanced (on opposite sides of mean)
			// balance them with Bernoulli probability equal to 1 - imbalance parameter
			else {
				if (!bdist(generateRandom())) {
					competition_fecundity[i][j] = fecundity_mean - fecundity_fecundity_diff_upper;
				}
			}
			if (growth_fecundity_diff_upper * growth_fecundity_diff_lower > 0.) {
				if (bdist(generateRandom())) {
					competition_growth[i][j] = growth_mean - growth_fecundity_diff_upper;
				}
			}
			else {
				if (!bdist(generateRandom())) {
					competition_growth[i][j] = growth_mean - growth_fecundity_diff_upper;
				}
			}
		}
	}
	return;
}


void Simulation::getFecundityGrowthCorrelation() {
	/* calculate Pearson's product moment correlation of interactions between growth and fecundity competition matrices.
	the correlation of these matrices is a metric of whether there are trade-offs between growth and reproduction. 
	if positive, species good at growing are also good at reproducing, and if negative, there is a trade-off. */
	
	int i, j;
	double covariance = 0.;
	double fecundity_variance = 0.;
	double growth_variance = 0.;
	double fecundity_mean = 0.;
	double growth_mean = 0.;

	// calculate mean of fecundity and growth competition matrices
	for (i = 0; i < num_species; i++) {
		for (j = 0; j < num_species; j++) {
			fecundity_mean += competition_fecundity[i][j];
			growth_mean += competition_growth[i][j];
		}
	}
	fecundity_mean = fecundity_mean / num_species / num_species;
	growth_mean = growth_mean / num_species / num_species;

	// calculate covariate and variance of fecundity and growth competition matrices
	for (i = 0; i < num_species; i++) {
		for (j = 0; j < num_species; j++) {
			covariance += (competition_fecundity[i][j] - fecundity_mean) * (competition_growth[i][j] - growth_mean);
			fecundity_variance += pow(competition_fecundity[i][j] - fecundity_mean, 2);
			growth_variance += pow(competition_growth[i][j] - growth_mean, 2);
		}
	}
	covariance = covariance / (num_species * num_species - 1.);
	fecundity_variance = fecundity_variance / (num_species * num_species - 1.);
	growth_variance = growth_variance / (num_species * num_species - 1.);

	// calculate correlation by dividing by standard deviation
	fecundity_growth_correlation = covariance / sqrt(fecundity_variance * growth_variance);

	return;
}


void Simulation::getImbalanceMean() {
	/* calculate imbalance of paired elements in growth an fecundity competition matrices.
	will return one value of mean imbalance for each competition type. imbalance occurs between
	species pairs when the first species has a larger effect on the second than the second has on 
	the first, relative to the mean strength of interactions in the community. */

	int i, j;
	fecundity_imbalance_mean = 0.;
	growth_imbalance_mean = 0.;

	// calculate mean of the differences between paired elements (imbalance)
	for (i = 0; i < num_species; i++) {
		for (j = i + 1; j < num_species; j++) {
			if (i != j) {
				fecundity_imbalance_mean += fabs(competition_fecundity[i][j] - competition_fecundity[j][i]);
				growth_imbalance_mean += fabs(competition_growth[i][j] - competition_growth[j][i]);
			}
		}
	}
	fecundity_imbalance_mean = fecundity_imbalance_mean / num_species / (num_species - 1.);
	growth_imbalance_mean = growth_imbalance_mean / num_species / (num_species - 1.);

	return;
}

void Simulation::getDiscreteTransitivity() {
	/* calculate discrete transitivity of growth and fecundity competition matrices.
	discrete transitivity is measured as  the relative intransitivity index, a metric
	that uses the competitive outcomes matrix (a binary matrix). completely transitive matrices 
	have RI index of 0 and maximally intransitive matrices have RI index of 1. */

	int i, j;
	double top_mean = 0.;
	double bottom_mean = 0.;
	double top_variance = 0.;
	double bottom_variance = 0.;
	double fecundity_rank_mean = 0.;	
	double growth_rank_mean = 0.;
	double fecundity_variance = 0.;
	double growth_variance = 0.;

	int *top_rank = new int[num_species];
	int *bottom_rank = new int[num_species];
	if (!top_rank || !bottom_rank ) {
		fprintf(stderr, "Error, unable to allocate memory for transitivity calculation\n");
		exit(-1);
	}

	for (i = 0; i < num_species; i++ ) {
		fecundity_rank[i] = 0.;
		growth_rank[i] = 0.;
	}

	// for a community with the same  number of species, top_rank describes the number of species
	// that each species outcompetes in a completely transitive community.
	// bottom_rank describes the number of species that each species outcompetes in a 
	// maximally intransitive community. the variance of ranks for these two cases
	// represents the maximum and minimum variance of ranks used in the RI index
	for (i = 0; i < num_species; i++) {
		top_rank[i] = num_species - (i + 1);
		if (num_species % 2 == 0) {
			if (i < num_species / 2)
				bottom_rank[i] = num_species / 2;
			else
				bottom_rank[i] = num_species / 2. - 1;
		}
		else
			bottom_rank[i] = (num_species - 1) / 2;

		top_mean += top_rank[i];
		bottom_mean += bottom_rank[i];
	}
	top_mean = top_mean / num_species;
	bottom_mean = bottom_mean / num_species;

	for (i = 0; i < num_species; i++) {
		top_variance += pow(top_rank[i] - top_mean, 2);
		bottom_variance += pow(bottom_rank[i] - bottom_mean, 2);
	}
	top_variance = top_variance / (num_species - 1);
	bottom_variance = bottom_variance / (num_species - 1);

	delete[] top_rank;
	delete[] bottom_rank;
	
	// create competitive outcomes matrix for growth and fecundity competition matrices. 
	// for matrix index i, j: 1 if j outcompetes i, and 0 if i outcompetes j
	for (i = 0; i < num_species; i++) {
		for (j = i + 1; j < num_species; j++) {

			if (competition_fecundity[i][j] < competition_fecundity[j][i]) {
				fecundity_transitivity[i][j] = 1.;
				fecundity_transitivity[j][i] = 0.;
			}
			else if (competition_fecundity[i][j] > competition_fecundity[j][i]) {
				fecundity_transitivity[i][j] = 0.;
				fecundity_transitivity[j][i] = 1.;
			}
			fecundity_rank[i] += fecundity_transitivity[j][i];
			fecundity_rank[j] += fecundity_transitivity[i][j];

			if (competition_growth[i][j] < competition_growth[j][i]) {
				growth_transitivity[i][j] = 1.;
				growth_transitivity[j][i] = 0.;
			}
			else if (competition_growth[i][j] > competition_growth[j][i]) {
				growth_transitivity[i][j] = 0.;
				growth_transitivity[j][i] = 1.;
			}
			growth_rank[i] += growth_transitivity[j][i];
			growth_rank[j] += growth_transitivity[i][j];
		}	
	}

	// calculate variance of ranks for the observed community
	// and compare to minimum and maximum variance in RI index calculation.
	for (i = 0; i < num_species; i++) {
		fecundity_rank_mean += fecundity_rank[i];
		growth_rank_mean += growth_rank[i];
	}
	fecundity_rank_mean = fecundity_rank_mean / num_species;
	growth_rank_mean = growth_rank_mean / num_species;

	for (i = 0; i < num_species; i++) {
		fecundity_variance += pow(fecundity_rank[i] - fecundity_rank_mean, 2);
		growth_variance += pow(growth_rank[i] - growth_rank_mean, 2);
	}
	fecundity_variance = fecundity_variance / (num_species - 1);
	growth_variance = growth_variance / (num_species - 1);

	fecundity_relative_intransitivity = 1. - (fecundity_variance - bottom_variance) / (top_variance - bottom_variance);
	growth_relative_intransitivity = 1. - (growth_variance - bottom_variance) / (top_variance - bottom_variance);

	return;
}


void Simulation::getContinuousTransitivity() {
	/* calculate continous transitivity of growth and fecundity competition matrices.
	continous transitivity uses the sum of pairwise elements scaled by the range over which
	competition occurs. the magnitude is equal for pairwise elements, but the species that 
	is a better competitor is positive and the worse competitor is negative. ranks are calculated
	as with discrete transitivity. this metric is not currently implemented. */

	int i, j;

	for (i = 0; i < num_species; i++) {
		for (j = i + 1; j < num_species; j++) {
			fecundity_transitivity[i][j] = (-competition_fecundity[i][j] + competition_fecundity[j][i])
											/ (competition_upper_bound - competition_lower_bound);
			fecundity_transitivity[j][i] = -fecundity_transitivity[i][j];
			fecundity_rank[i] += fecundity_transitivity[j][i];
			fecundity_rank[j] += fecundity_transitivity[i][j];

			growth_transitivity[i][j] = (-competition_growth[i][j] + competition_growth[j][i])
											/ (competition_upper_bound - competition_lower_bound);

			growth_transitivity[j][i] = -growth_transitivity[i][j];
			growth_rank[i] += growth_transitivity[j][i];
			growth_rank[j] += growth_transitivity[i][j];
		}	
	}

	return;
}


void Simulation::getDiscreteFecundityTransitivity() {
	/* calculate the ranks of the fecundity competition matrix. used when manipulating transitivity
	in method 'setCompetitionTransitivity'. */

	int i, j, k;
	double pivot = 0;

	// create competitive outcomes matrix (fecundity_transitivity), for matrix index i, j: 1 if j
	// outcompetes i, and 0 if i outcompetes j. calculate ranks, or the number of species that each species outcompetes.
	for (i = 0; i < num_species; i++) {
		 fecundity_rank[i] = 0.;
	}

	for (i = 0; i < num_species; i++) {
		for (j = i + 1; j < num_species; j++) {

			if (competition_fecundity[i][j] < competition_fecundity[j][i]) {
				fecundity_transitivity[i][j] = 1.;
				fecundity_transitivity[j][i] = 0.;
			}
			else if (competition_fecundity[i][j] > competition_fecundity[j][i]) {
				fecundity_transitivity[i][j] = 0.;
				fecundity_transitivity[j][i] = 1.;
			}
			fecundity_rank[i] += fecundity_transitivity[j][i];
			fecundity_rank[j] += fecundity_transitivity[i][j];
		}	
	}
	
	fprintf(stdout, "fecundity_rank: ");
	for (i = 0; i < num_species; i++) {
		fprintf(stdout, "%.0f ", fecundity_rank[i]);
	}
	fprintf(stdout, "\n");

	return;
}


void Simulation::getDiscreteGrowthTransitivity() {
	/* calculate the ranks of the growth competition matrix. used when manipulating transitivity
	in method 'setCompetitionTransitivity'. */

	int i, j, k;
	double pivot = 0;

	// create competitive outcomes matrix (fecundity_transitivity), for matrix index i, j: 1 if j
	// outcompetes i, and 0 if i outcompetes j. calculate ranks, or the number of species that each species outcompetes.
	for (i = 0; i < num_species; i++) {
		 growth_rank[i] = 0.;
	}
	for (i = 0; i < num_species; i++) {
		for (j = i + 1; j < num_species; j++) {
			if (competition_growth[i][j] < competition_growth[j][i]) {
				growth_transitivity[i][j] = 1.;
				growth_transitivity[j][i] = 0.;
			}
			else if (competition_growth[i][j] > competition_growth[j][i]) {
				growth_transitivity[i][j] = 0.;
				growth_transitivity[j][i] = 1.;
			}
			growth_rank[i] += growth_transitivity[j][i];
			growth_rank[j] += growth_transitivity[i][j];
		}	
	}
	
	fprintf(stdout, "growth_rank: ");
	for (i = 0; i < num_species; i++) {
		fprintf(stdout, "%.0f ", growth_rank[i]);
	}
	fprintf(stdout, "\n");

	return;
}


void Simulation::setCompetitionTransitivity() {
	/* manipulate the transitivity of fecundity and growth competition matrices. this is controlled by the
	parameters 'fecundity_transitivity_type' and 'growth_transitivity_type.' if these are set to 1, 
	the competition matrices will be completely transitive. if both are set to -1, the competition matrices
	will be maximally intransitive. may call on method 'setGrowthCompetitionTransitivity' if transitivity
	types are not the same */

	int i, j, k;
	int ih, jh, complete;
	double hold;
	int *desired_fecundity_rank = new int[num_species];
	if(!desired_fecundity_rank) {
		fprintf(stderr, "Error, unable to allocate memory to initialize transitive competition\n");
		exit(-1);
	}

	// set desired pattern of ranks (number of species outcompeted)
	// completely transitive
	if (fecundity_transitivity_type == 1.) {
		fprintf(stdout, "fecundity_transitivity_type == 1. (transitive)\n");
		for (i = 0; i < num_species; i++) {
			desired_fecundity_rank[i] = num_species - (i + 1);
		}		
	}
	// maximally intransitive
	else if (fecundity_transitivity_type == -1.) {
		fprintf(stdout, "fecundity_transitivity_type == -1. (intransitive)\n");

		if (num_species % 2 == 0) {
			for (i = 0; i < num_species / 2; i++) {
				desired_fecundity_rank[i] = num_species / 2;
				desired_fecundity_rank[num_species - (i + 1)] = num_species / 2 -1; 
			}
		}
		else {
			for (i = 0; i < num_species; i++) {
				desired_fecundity_rank[i] = (num_species - 1) / 2;
			}
		}
	}

	// if fecundity transitivity type is set to be random, but growth transitivity type is not, 
	// call method to set growth transitivity
	else if (fecundity_transitivity_type == 0) {
		getDiscreteFecundityTransitivity();
		setGrowthCompetitionTransitivity();
		return;
	}

	fprintf(stdout, "desired_fecundity_rank: ");
	for (i = 0; i < num_species; i++) {
		fprintf(stdout, "%d ", desired_fecundity_rank[i]);
	}
	fprintf(stdout, "\n");

  	// calculate ranks of observed fecundity matrix, to be compared to desired ranks
	getDiscreteFecundityTransitivity();

	int * fecundity_diff = new int[num_species];
	if (!fecundity_diff) {
		fprintf(stderr, "Error, unable to allocate memory to set transitivity difference array\n");
		exit(-1);
	}

	// calculate the difference between the observed ranks and the desired ranks
	// fecundity_diff[i]: number of additional species that species i needs to beat if negative
	// number of additional species that species i needs to lose to if positive
	for (i = 0; i < num_species; i++) {
		fecundity_diff[i] = (int) desired_fecundity_rank[i] - fecundity_rank[i];
	}

	fprintf(stdout, "fecundity_diff: ");
	for (i = 0; i < num_species; i++)
		fprintf(stdout, "%d ", fecundity_diff[i]);
	fprintf(stdout, "\n");

	// switch pairwise effects until observed ranks equal desired values
	for (i = 0; i < num_species; i++) {
		for (j = num_species - 1; j > i; j--) {
			// if species i needs to outcompete more species and species j needs to compete fewer species, switch pairwise effects
			if (fecundity_diff[i] < 0 && fecundity_transitivity[i][j] == 0) {
				hold = competition_fecundity[j][i];
				competition_fecundity[j][i] = competition_fecundity[i][j];
				competition_fecundity[i][j] = hold;
				fecundity_diff[i]++;
				fecundity_diff[j]--;
			}
			// if species i needs to outcompete fewer species and species j needs to compete more species, switch pairwise effects
			else if (fecundity_diff[i] > 0 && fecundity_transitivity[i][j] == 1) {
				hold = competition_fecundity[j][i];
				competition_fecundity[j][i] = competition_fecundity[i][j];
				competition_fecundity[i][j] = hold;
				fecundity_diff[i]--;
				fecundity_diff[j]++;
			}
		}
	}

	// in intransitive case, a single transitive triplet could remain after this process, check for a difference between desired and observed ranks
	ih = -1;
	jh = -1;
	while (1) {
		complete = 0;
		getDiscreteFecundityTransitivity();  // calculate ranks of manipulated fecundity competition matrix

		for (i = 0; i < num_species; i++) {
			fecundity_diff[i] = (int) desired_fecundity_rank[i] - fecundity_rank[i];
		}

		fprintf(stdout, "fecundity_diff: ");
		for (i = 0; i < num_species; i++)
			fprintf(stdout, "%d ", fecundity_diff[i]);
		fprintf(stdout, "\n");
		fprintf(stdout, "ih = %d, jh = %d\n", ih, jh);

		for (i = 0; i < num_species; i++) {
			if (fecundity_diff[i] > 0)
				ih = i;
			else if (fecundity_diff[i] < 0)
				jh = i;
			else if (fecundity_diff[i] == 0)
				complete++;

		}	

		if (complete == num_species)
			break;

		// if transitive triplet remains, switch third species and allow the unresolved pair to flip
		for (i = 0; i < num_species; i++) {
			if (i == ih || i == jh)
				continue;
			if (fecundity_transitivity[ih][i] == 1 && fecundity_transitivity[jh][i] == 0) {
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
	
	getDiscreteFecundityTransitivity();
	getDiscreteGrowthTransitivity();

	delete[] desired_fecundity_rank;
	delete[] fecundity_diff;

	// calculate ranks of manipulated fecundity competition matrix
	getDiscreteTransitivity();

	// if growth competition matrix has the same transitivity type as the fecundity matrix,	use the fecundity matrix (already set)
	// to manipulate elements in the growth matrix. if both growth and fecundity competition matrices are intransitive, 
	// they can have either the same species hierarchy or they can have reverse hierarchies (growth/fecundity trade-off). 
	if (growth_transitivity_type == fecundity_transitivity_type) {
		for (i = 0; i < num_species; i++) {
			for (j = i + 1; j < num_species; j++) {
				if (i == j)
					continue;
				if (fecundity_transitivity[i][j] != growth_transitivity[i][j] && fecundity_growth_relative_hierarchy != -1) {
					hold = competition_growth[j][i];
					competition_growth[j][i] = competition_growth[i][j];
					competition_growth[i][j] = hold;
				}
				else if (fecundity_transitivity[i][j] == growth_transitivity[i][j] && fecundity_growth_relative_hierarchy == -1) {
					hold = competition_growth[j][i];
					competition_growth[j][i] = competition_growth[i][j];
					competition_growth[i][j] = hold;
				}
			}
		}
		if (fecundity_growth_relative_hierarchy == 0) {
			shuffleMatrix(competition_growth, num_species);
		}
	}
	else if (growth_transitivity_type != 0) {
		setGrowthCompetitionTransitivity();
	}

	getDiscreteGrowthTransitivity();

	return;
}

void Simulation::shuffleArray(int *array, int n) {
    if (n > 1) {
        int i;
        for (i = 0; i < n - 1; i++) {
          int j = i + getRandom() / (global_random_generator.max() / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
	return;
}

void Simulation::shuffleMatrix(double **array, int n) {
	int i, j;

	double ** hold_matrix = new double*[n];
	for (i = 0; i < n; i++) {
		hold_matrix[i] = new double[n];
		for (j = 0; j < n; j++) {
			hold_matrix[i][j] = array[i][j];
		}
	}

	int * shuffle_index = new int[n];
	for (i = 0; i < n; i++) {
		shuffle_index[i] = i;
	}

	shuffleArray(shuffle_index, n);

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			array[shuffle_index[i]][shuffle_index[j]] = hold_matrix[i][j];
		}
	}	

	for (i = 0; i < n; i++) {
		delete [] hold_matrix[i];
	}
	delete [] hold_matrix;

	delete [] shuffle_index;

	return;
}


void Simulation::setGrowthCompetitionTransitivity() {
	/* this method is called in 'setCompetitionTransitivity' after the fecundity transitivity matrix is set, 
	and only if the growth transitivity type is not the same as fecundity. this method sets growth competition matrix
	 transitivity according to parameter 'growth_transitivity_type.' */

	int i, j, k;
	int ih, jh, complete;
	double hold;
	int *desired_growth_rank = new int[num_species];
	if (!desired_growth_rank) {
		fprintf(stderr, "Error, unable to allocate memory to initialize transitive competition\n");
		exit(-1);
	}
	// set desired pattern of ranks
	// completely transitive
	if (growth_transitivity_type == 1.) {
		fprintf(stdout, "growth_transitivity_type == 1. (transitive)\n");
		for (i = 0; i < num_species; i++)
			desired_growth_rank[i] = num_species - (i + 1);
	}
	// maximally intransitive
	else if (growth_transitivity_type == -1.) {
		fprintf(stdout, "growth_transitivity_type == -1. (intransitive)\n");
		if (num_species % 2 == 0) {
			for (i = 0; i < num_species / 2; i++) {
				desired_growth_rank[i] = num_species / 2 ;
				desired_growth_rank[num_species - (i + 1)] = num_species / 2 - 1 ; 
			}
		}
		else {
			for (i = 0; i < num_species; i++) {
				desired_growth_rank[i] = (num_species - 1) / 2 ;
			}
		}
	}

	fprintf(stdout, "desired_growth_rank: ");
	for (i = 0; i < num_species; i++) {
		fprintf(stdout, "%d ", desired_growth_rank[i]);
	}
	fprintf(stdout, "\n");

  	// calculate ranks of observed growth matrix, to be compared to desired ranks
	getDiscreteGrowthTransitivity();

	int * growth_diff = new int[num_species];
	if (!growth_diff) {
		fprintf(stderr, "Error, unable to allocate memory for setting transitivity difference array\n");
		exit(-1);
	}

	// calculate the difference between the observed ranks and the desired ranks
	// growth_diff[i]: number of additional species that species i needs to beat if negative
	// number of additional species that species i needs to lose to if positive
	//
	for (i = 0; i < num_species; i++) {
		growth_diff[i] = (int) desired_growth_rank[i] - growth_rank[i];
	}
	

	fprintf(stdout, "growth_diff: ");
	for (i = 0; i < num_species; i++)
		fprintf(stdout, "%d ", growth_diff[i]);
	fprintf(stdout, "\n");

	// switch pairwise effects until observed ranks equal desired values
	for (i = 0; i < num_species; i++) {
		for (j = num_species - 1; j > i; j--) {
			// if species i needs to outcompete more species and species j needs to compete fewer species, switch pairwise effects
			if (growth_diff[i] < 0 && growth_transitivity[i][j] == 0) {
				hold = competition_growth[j][i];
				competition_growth[j][i] = competition_growth[i][j];
				competition_growth[i][j] = hold;
				growth_diff[i]++;
				growth_diff[j]--;
			}
			// if species i needs to outcompete fewer species and species j needs to compete more species, switch pairwise effects
			else if (growth_diff[i] > 0 && growth_transitivity[i][j] == 1) {
				hold = competition_growth[j][i];
				competition_growth[j][i] = competition_growth[i][j];
				competition_growth[i][j] = hold;
				growth_diff[i]--;
				growth_diff[j]++;
			}
		}
	}

	// in intransitive case, a single transitive triplet could remain after this process, check for a difference between desired and observed ranks
	ih = -1;
	jh = -1;
	while (1) {
		complete = 0;
		getDiscreteGrowthTransitivity();  // calculate ranks of manipulated growth competition matrix

		for (i = 0; i < num_species; i++) {
			growth_diff[i] = (int) desired_growth_rank[i] - growth_rank[i];
		}

		for (i = 0; i < num_species; i++) {
			if (growth_diff[i] > 0)
				ih = i;
			else if (growth_diff[i] < 0)
				jh = i;
		}	
		
		fprintf(stdout, "growth_diff: ");
		for (i = 0; i < num_species; i++)
			fprintf(stdout, "%d ", growth_diff[i]);
		fprintf(stdout, "\n");
		fprintf(stdout, "ih = %d, jh = %d\n", ih, jh);
		
		for (i = 0; i < num_species; i ++) {
			if (growth_diff[i] == 0)
				complete++;
		}
		if (complete == num_species)
			break;

		// if transitive triplet remains, switch third species and allow the unresolved pair to flip
		for (i = 0; i < num_species; i++) {
			if (i == ih || i == jh)
				continue;
			if (growth_transitivity[ih][i] == 1 && growth_transitivity[jh][i] == 0) {
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

	shuffleMatrix(competition_growth, num_species);

	delete[] desired_growth_rank;
	delete[] growth_diff;

	return;
}


double Simulation::getRandomUniformReal(double lower_bound, double upper_bound) {

	return lower_bound + (upper_bound - lower_bound) * ( (double) getRandom() ) / ( (double) global_random_generator.max() );

}

double Simulation::getRandomNormal(double mean, double sdev) {

	static bool hasSpare = false;
	static double spare;
	double u, v; 
	double s = 0.;

	if(hasSpare) {
		hasSpare = false;
		return mean + sdev * spare;
	}

	hasSpare = true;
	while( (s >=1.0) || (s==0) ) {
		u = ( (double) getRandom() ) / ( (double) global_random_generator.max() ) * 2.0 - 1.0;
		v = ( (double) getRandom() ) / ( (double) global_random_generator.max() ) * 2.0 - 1.0;
		s = u * u + v * v;
	}

	s = sqrt(-2.0 * log(s) / s);
	spare = v * s;

	return mean + sdev * u * s;

}


