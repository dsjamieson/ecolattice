
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2018			
	 *
	 *		methods for the Simulation class. sets up growth and 
	 *		fecundity competition matrices; manipulates imbalance,
	 *		transitivity, and correlation as specified by parameters.
	 *
	 ***********************************************************/

#include "ecolattice.hpp"

void Ecolattice::initializeUniformCompetition(void){
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


void Ecolattice::initializeTNormalCompetition(void) {
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


void Ecolattice::initializeUniformCorrelatedCompetition(void){
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


void Ecolattice::initializeTNormalCorrelatedCompetition(void){
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

void Ecolattice::imbalanceCompetition(void) {
	/* manipulate imbalance of (already initialized) growth and fecundity competition matrices.
	imbalance is governed by 'imbalance' parameter, the Bernoulli probability of increasing
	the imbalance of paired elements. imbalance refers to one species having a greater effect 
	on a second species than the second has on the first. */

	int i, j;
	double fecundity_fecundity_diff_upper, fecundity_fecundity_diff_lower, growth_fecundity_diff_upper, growth_fecundity_diff_lower;
	double fecundity_mean = 0;
	double growth_mean = 0;
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

void Ecolattice::setCompetitionTransitivity(void) {
	/* manipulate the transitivity of fecundity and growth competition matrices. this is controlled by the
	parameters 'fecundity_transitivity_type' and 'growth_transitivity_type.' if these are set to 1, 
	the competition matrices will be completely transitive. if both are set to -1, the competition matrices
	will be maximally intransitive. may call on method 'setGrowthCompetitionTransitivity' if transitivity
	types are not the same */

	int i, j;
	int ih, jh, complete;
	double hold;
	std::vector<int> desired_fecundity_rank(num_species);

	// set desired pattern of ranks (number of species outcompeted)
	// completely transitive
	if (fecundity_transitivity_type == 1.) {
		//fprintf(stdout, "fecundity_transitivity_type == 1. (transitive)\n");
		for (i = 0; i < num_species; i++) {
			desired_fecundity_rank[i] = num_species - (i + 1);
		}		
	}
	// maximally intransitive
	else if (fecundity_transitivity_type == -1.) {
		//fprintf(stdout, "fecundity_transitivity_type == -1. (intransitive)\n");
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
	/*
	fprintf(stdout, "desired_fecundity_rank: ");
	for (i = 0; i < num_species; i++) {
		fprintf(stdout, "%d ", desired_fecundity_rank[i]);
	}
	fprintf(stdout, "\n");
	*/

  	// calculate ranks of observed fecundity matrix, to be compared to desired ranks
	getDiscreteFecundityTransitivity();

	std::vector<int>  fecundity_diff(num_species);

	// calculate the difference between the observed ranks and the desired ranks
	// fecundity_diff[i]: number of additional species that species i needs to beat if negative
	// number of additional species that species i needs to lose to if positive
	for (i = 0; i < num_species; i++) {
		fecundity_diff[i] = (int) desired_fecundity_rank[i] - fecundity_rank[i];
	}
	/*
	fprintf(stdout, "fecundity_diff: ");
	for (i = 0; i < num_species; i++)
		fprintf(stdout, "%d ", fecundity_diff[i]);
	fprintf(stdout, "\n");
	*/

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

		/*
		fprintf(stdout, "fecundity_diff: ");
		for (i = 0; i < num_species; i++)
			fprintf(stdout, "%d ", fecundity_diff[i]);
		fprintf(stdout, "\n");
		fprintf(stdout, "ih = %d, jh = %d\n", ih, jh);
		*/

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
			shuffleMatrix(competition_growth);
		}
	}
	else if (growth_transitivity_type != 0) {
		setGrowthCompetitionTransitivity();
	}

	getDiscreteGrowthTransitivity();

	return;
}


void Ecolattice::setGrowthCompetitionTransitivity(void) {
	/* this method is called in 'setCompetitionTransitivity' after the fecundity transitivity matrix is set, 
	and only if the growth transitivity type is not the same as fecundity. this method sets growth competition matrix
	 transitivity according to parameter 'growth_transitivity_type.' */

	int i, j;
	int ih, jh, complete;
	double hold;
	std::vector<int> desired_growth_rank(num_species);
	// set desired pattern of ranks
	// completely transitive
	if (growth_transitivity_type == 1.) {
		//fprintf(stdout, "growth_transitivity_type == 1. (transitive)\n");
		for (i = 0; i < num_species; i++)
			desired_growth_rank[i] = num_species - (i + 1);
	}
	// maximally intransitive
	else if (growth_transitivity_type == -1.) {
		//fprintf(stdout, "growth_transitivity_type == -1. (intransitive)\n");
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

	/*
	fprintf(stdout, "desired_growth_rank: ");
	for (i = 0; i < num_species; i++) {
		fprintf(stdout, "%d ", desired_growth_rank[i]);
	}
	fprintf(stdout, "\n");
	*/

  	// calculate ranks of observed growth matrix, to be compared to desired ranks
	getDiscreteGrowthTransitivity();

	std::vector<int> growth_diff(num_species);

	// calculate the difference between the observed ranks and the desired ranks
	// growth_diff[i]: number of additional species that species i needs to beat if negative
	// number of additional species that species i needs to lose to if positive
	//
	for (i = 0; i < num_species; i++) {
		growth_diff[i] = (int) desired_growth_rank[i] - growth_rank[i];
	}
	
	/*
	fprintf(stdout, "growth_diff: ");
	for (i = 0; i < num_species; i++)
		fprintf(stdout, "%d ", growth_diff[i]);
	fprintf(stdout, "\n");
	*/

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
		/*
		fprintf(stdout, "growth_diff: ");
		for (i = 0; i < num_species; i++)
			fprintf(stdout, "%d ", growth_diff[i]);
		fprintf(stdout, "\n");
		fprintf(stdout, "ih = %d, jh = %d\n", ih, jh);
		*/
		
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
	shuffleMatrix(competition_growth);
	return;
}

void Ecolattice::shuffleArray(std::vector<int> & t_vector) {
	/* randomly shuffle elements of a 1D vector. used when setting transitivity  */

    if (t_vector.size() > 1) {
        unsigned long i;
        for (i = 0; i < t_vector.size() - 1; i++) {
          int j = i + getRandom() / (static_cast<double>(std::mt19937::max()) / (t_vector.size() - i) + 1);
          int t = t_vector[j];
          t_vector[j] = t_vector[i];
          t_vector[i] = t;
        }
    }
	return;
}

void Ecolattice::shuffleMatrix(std::vector<std::vector<double>> & t_vector) {
	/* randomly shuffle elements of a 2D vector. used when setting transitivity. */

	unsigned long i, j;

	std::vector<std::vector<double>> hold_matrix(t_vector.size());
	for (i = 0; i < t_vector.size(); i++) {
		hold_matrix[i].resize(t_vector[i].size());
		for (j = 0; j < t_vector.size(); j++) {
			hold_matrix[i][j] = t_vector[i][j];
		}
	}

	std::vector<int> shuffle_index(t_vector.size());
	for (i = 0; i < t_vector.size(); i++) {
		shuffle_index[i] = i;
	}

	shuffleArray(shuffle_index);

	for (i = 0; i < t_vector.size(); i++) {
		for (j = 0; j < t_vector.size(); j++) {
			t_vector[shuffle_index[i]][shuffle_index[j]] = hold_matrix[i][j];
		}
	}	
	return;
}
