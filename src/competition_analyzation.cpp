
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2020			
	 *
	 *		methods for the Ecolattice class: these methods do not 
	 *		set any values. calculates competition characteristics:
	 *		abundance, imbalance, transitivity, and correlation.
	 *
	 ***********************************************************/

#include "ecolattice.hpp"

void Ecolattice::getSpeciesAbundance(void) {
	/* calculate the number of individuals of each species in the lattice. this only happens after initialization, before
  	   any simulation time steps. subsequently, species abundances are monitored using the incrementSpeciesAbundance
  	   and decrementSpeciesAbundance methods in the SiteStepper class. */
	for (int i = 0; i < num_species; i++)
		species_abundance[i] = 0;
	#pragma omp parallel for
	for (int i = 0; i < lattice_size; i++) {
		for (int j = 0; j < lattice_size; j++) {
			int s = abs(lattice[i][j]);
			if (s != 0) {
				#pragma omp atomic
				species_abundance[s - 1]++;
			}
		}
	}
	return;
}

void Ecolattice::getFecundityGrowthCorrelation(void) {
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


void Ecolattice::getImbalanceMean(void) {
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

void Ecolattice::getDiscreteTransitivity(void) {
	/* calculate transitivity of growth and fecundity competition matrices.
	transitivity is measured with the relative intransitivity index, a metric
	that uses the competitive outcomes matrix (a binary matrix). completely transitive matrices 
	have an RI index of 0 and maximally intransitive matrices have RI index of 1. */

	int i, j;
	double top_mean = 0.;
	double bottom_mean = 0.;
	double top_variance = 0.;
	double bottom_variance = 0.;
	double fecundity_rank_mean = 0.;	
	double growth_rank_mean = 0.;
	double fecundity_variance = 0.;
	double growth_variance = 0.;

	std::vector<int> top_rank(num_species);
	std::vector<int> bottom_rank(num_species);
	for (i = 0; i < num_species; i++ ) {
		fecundity_rank[i] = 0.;
		growth_rank[i] = 0.;
	}

	// for a community with the same number of species, top_rank describes the number of species
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

void Ecolattice::getDiscreteFecundityTransitivity(void) {
	/* calculate the ranks of the fecundity competition matrix. used when manipulating transitivity
	in method 'setCompetitionTransitivity'. */

	int i, j;

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
	/*
	fprintf(stdout, "fecundity_rank: ");
	for (i = 0; i < num_species; i++) {
		fprintf(stdout, "%.0f ", fecundity_rank[i]);
	}
	fprintf(stdout, "\n");
	*/

	return;
}


void Ecolattice::getDiscreteGrowthTransitivity(void) {
	/* calculate the ranks of the growth competition matrix. used when manipulating transitivity
	in method 'setCompetitionTransitivity'. */

	int i, j;

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
	/*
	fprintf(stdout, "growth_rank: ");
	for (i = 0; i < num_species; i++) {
		fprintf(stdout, "%.0f ", growth_rank[i]);
	}
	fprintf(stdout, "\n");
	*/
	return;
}

void Ecolattice::getContinuousTransitivity(void) {
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
