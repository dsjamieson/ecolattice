
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2018			
	 *
	 *		methods for the Simulation class, for omp version
	 *		of the code.
	 *
	 ***********************************************************/

#include "simulation.h"


void Simulation::seedGenerator(std::mt19937& this_random_generator) {
	/* create a random vector of seeds (a seed sequence) given the seeds specified randomly or in
	the parameter file (for restarted simulations). seeds fed to the global RNG. */
    std::seed_seq seq(seeds, seeds + 5);
	std::vector<std::uint32_t> seed_vector(std::mt19937::state_size);
    seq.generate(seed_vector.begin(), seed_vector.end());
	std::seed_seq seq2(seed_vector.begin(), seed_vector.end());
	this_random_generator.seed(seq2);
	return;
}

std::mt19937& Simulation::generateRandom(unsigned long long& this_random_count, std::mt19937& this_random_generator) {
	/* add to the random count and get a random draw from the global RNG. */

	this_random_count += 2;
	return this_random_generator;
}

void Simulation::discardRandom(unsigned long long n, unsigned long long& this_random_count, std::mt19937& this_random_generator) {
	/* add to the random count and discard values from the global RNG. */

	this_random_generator.discard(n);
	this_random_count += n;
}

unsigned long long Simulation::getMaxRandomCount() {
	return max_random_count;
}



void Simulation::updateSingleSite(int i, int j, unsigned long long& this_random_count, std::mt19937& this_random_generator) {
	/* this method runs through all processes, including germination, survival, growth, reproduction, and death,
	for a single site in the lattice. workers use this method to update their local copies of 'next_lattice' and 'next_dispersal_lattice' */
	fprintf(stdout, "site = %d, %d\n", i, j);
	int k, l, kp, lp;
	unsigned long long start_random_count = this_random_count;

	int this_delta = 0;
	int this_neighborhood_size = 0;
	int this_species = 0;
	int this_stage = 0;
	int this_max_dispersal = 0;
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

	this_species = abs(lattice[i][j]);  // species in this site i, j
	fprintf(stdout, "this_species = %d\n", this_species);

	// if the site is open, determine whether or not something germinates with a Bernoulli probability based on the total number of seeds in the site
	// if something germinates, select the species that germinates with probability equal to the relative abundance of each species's seeds in this site.
	if (this_species == 0) {
	
		for (k = 0; k < num_species; k++)
			total_seeds += dispersal_lattice[i][j][k];
		if (i == 9 && j == 7) {
			fprintf(stdout, "total_seeds = %f\n", total_seeds);
		}
		std::bernoulli_distribution germ_dist(1. - pow(1. - germination_probability, total_seeds));

		if (germ_dist(generateRandom(this_random_count, this_random_generator))) {
			std::discrete_distribution<int> species_dist(dispersal_lattice[i][j], dispersal_lattice[i][j] + num_species);
			l = abs(species_dist(generateRandom(this_random_count, this_random_generator)));
			next_lattice[i][j] = -(l + 1);
			if (i == 9 && j == 7) {
				fprintf(stdout, "next_lattice[9, 7] = %d\n", next_lattice[i][j]);
			}

			#pragma omp atomic
			species_abundance[l]++;
		}
		else {
			next_lattice[i][j] = 0;
		}
	}
	else {
		// if the site is not open, determine whether it's occupied by a juvenile or an adult
		this_stage = lattice[i][j] / abs(lattice[i][j]);
		// the individual survives with stage-specific probability
		if (this_stage > 0)
			this_survival_probability = adult_survival_probability[this_species - 1] ;
		else
			this_survival_probability = juvenile_survival_probability[this_species - 1];

		std::bernoulli_distribution survival_dist(this_survival_probability);
		// if species survives, it persists to next time step
		if (survival_dist(generateRandom(this_random_count, this_random_generator))) {
			next_lattice[i][j] = lattice[i][j];

			// determine abundance of each species in the neighborhood of the focal individual
			// neighborhood limits depend on the parameter delta for this_species
			this_delta = delta[this_species - 1];
			this_neighborhood_size = neighborhood_size[this_species - 1];
			// determine intrinsic fecundity, dispersal probability, dispersal length, and maximum competition based on the species identity
			this_intrinsic_fecundity = intrinsic_fecundity[this_species - 1];
			this_dispersal_probability = dispersal_probability[this_species - 1];
			this_dispersal_length = dispersal_length[this_species - 1];
			this_maximum_competition = maximum_competition[this_species - 1];

			neighborhood = new int[this_neighborhood_size + 1]; 
			if (!neighborhood) {
					fprintf(stderr, "Error, neighborhood memory allocation failed for site %d %d\n", i, j);	
					exit(-1);
			}
			for (k = 0; k < this_neighborhood_size + 1; k++)
				neighborhood[k] = 0;
			neighborhood_abundance = new int[num_species];
			if (!neighborhood_abundance) {
					fprintf(stderr, "Error, neighborhood abundance memory allocation failed for site %d %d\n", i, j);	
					exit(-1);
			}
			for (k = 0; k < num_species; k++)
				neighborhood_abundance[k] = 0;
			int neighbor_index = 0;	
			for (k = 0; k < 2 * this_delta + 1; k++) {
				for (l = 0; l < 2 * this_delta + 1; l++) {
					// skip center case
					if (k != this_delta || l != this_delta) {
						// lattice indices for site neighbors, accounting for periodic boundary conditions
						kp = k;
						lp = l;
						// if lattice row index is negative, wrap around to end of row
						if (i - this_delta + kp < 0)
							kp += lattice_size - i;
						// if lattice row index is past the end of the row, wrap around to start of row
						else if (i - this_delta + kp > lattice_size - 1)
							kp += -lattice_size;
						// if lattice column index is negative, wrap around to end of column
						if (j - this_delta + lp < 0)
							lp += lattice_size - j;
						// if lattice column index is past the end of the column, wrap around to start of column
						else if (j - this_delta + lp > lattice_size - 1)
							lp += -lattice_size;
						neighbor_index = l + (2 * this_delta + 1) * k;
						// use neighbor lattice indices to place neighbors in neighborhood vector
						neighborhood[neighbor_index] = lattice[i - this_delta + kp][j - this_delta + lp];
						if (neighborhood[neighbor_index] != 0) {
							neighborhood_abundance[abs(neighborhood[neighbor_index]) - 1]++;
							total_abundance++;
						}	
					}
				}
			}
			delete[] neighborhood;

			// if focal individual is a juvenile, it will grow to become an adult with a probability dependent on competition with individuals in its neighborhood
			// depends on the growth competition matrix, which dictates these interactions. the maximum probability is set as a parameter (never a 100% chance of growing)
			if (this_stage < 0) {
				for (k = 0; k < num_species; k++)
					growth_probability += competition_growth[this_species - 1][k] * ((double) neighborhood_abundance[k]) / ((double) this_neighborhood_size);
				growth_probability = exp(growth_probability);
				if (growth_probability > this_maximum_competition)
					growth_probability = this_maximum_competition;

				std::bernoulli_distribution stage_dist(growth_probability);

				if (stage_dist(generateRandom(this_random_count, this_random_generator))) {
					next_lattice[i][j] = abs(next_lattice[i][j]);
				}
			}
			else {
				// if focal individual is an adult, it will reproduce with a fecundity based on competition with individuals in its neighborhood
				if (total_abundance != 0) {
					for (k = 0; k < num_species; k++)
						this_fecundity += competition_fecundity[this_species - 1][k] * ((double) neighborhood_abundance[k]) / ((double) this_neighborhood_size);
					this_fecundity = this_intrinsic_fecundity * exp(this_fecundity);
				}
				else {
					this_fecundity = this_intrinsic_fecundity;
				}

				// fecundity, the number of seeds produced by the focal individuals, is spread over the entire matrix
				// the number of seeds at each site is an exponential function of Euclidean distance
				// probability of dispersal quickly decays with distance depending on the dispersal length parameter
				// dispersal can occur within two times the dispersal length of the species
				// setting a dispersal limit reduces size of distance_probability array and improves speed			
				this_max_dispersal = (int) round(2 * dispersal_length[this_species - 1]);
				distance_probability = new double*[this_max_dispersal + 1];
				if (!distance_probability) {
					fprintf(stderr, "Error, distance probability memory allocation failed for site %d %d\n", i, j);	
					exit(-1);
				}
				for (k = 0; k <= this_max_dispersal; k++) {
					distance_probability[k] = new double[this_max_dispersal];
					if (!distance_probability[k]) {
						fprintf(stderr, "Error, distance probability memory allocation failed for site %d %d\n", i, j);	
						exit(-1);
					}
					for (l = 1; l <= this_max_dispersal; l++) {
						double distance = sqrt((double) k * k + l * l);
						double probability = 0.;
						if (round(distance) <= this_max_dispersal) {
							probability = exp(log(this_dispersal_probability) / dispersal_length[this_species - 1] * distance);
							// every distance has four locations on lattice
							distance_probability_sum +=	4. * probability;
						}				
						distance_probability[k][l - 1] = probability;
					}
				}
				// start at lower right rectangular quadrant
				// determine the four lattice indices that correspond with each dispersal distance
				for (k = 0; k <= this_max_dispersal; k++) {
					for (l = 1; l <= this_max_dispersal; l++) {
						if (k == 0) {
							// in same row as site i, j (l = distance from center)
							// the four lattice indices are in the cardinal directions
							int i1 = (i + l) % lattice_size;
							int i2 = i - l;
							if (i2 < 0)
								i2 += lattice_size;
							int j1 = (j + l) % lattice_size;
							int j2 = j - l;
							if (j2 < 0)
								j2 += lattice_size;
							// same row, i1 to the right
							#pragma omp atomic
							next_dispersal_lattice[i][j1][this_species - 1] += this_fecundity * distance_probability[k][l - 1] / distance_probability_sum;
							// same row, i2 to the left
							#pragma omp atomic
							next_dispersal_lattice[i][j2][this_species - 1] += this_fecundity * distance_probability[k][l - 1] / distance_probability_sum;
							// same column, i1 up
							#pragma omp atomic
							next_dispersal_lattice[i1][j][this_species - 1] += this_fecundity * distance_probability[k][l - 1] / distance_probability_sum;
							// same column, i2 down
							#pragma omp atomic
							next_dispersal_lattice[i2][j][this_species - 1] += this_fecundity * distance_probability[k][l - 1] / distance_probability_sum;	
							}
						else {
							// square quadrants not in the same row as site i, j
							int i1 = (i + k) % lattice_size;
							int i2 = i - k;
							if (i2 < 0)
								i2 += lattice_size;
							int j1 = (j + l) % lattice_size;
							int j2 = j - l;
							if (j2 < 0)
								j2 += lattice_size;
							// upper right quadrant
							#pragma omp atomic
							next_dispersal_lattice[i1][j1][this_species - 1] += this_fecundity * distance_probability[k][l - 1] / distance_probability_sum;
							// upper left quadrant
							#pragma omp atomic
							next_dispersal_lattice[i2][j1][this_species - 1] += this_fecundity * distance_probability[k][l - 1] / distance_probability_sum;
							// lower right quadrant
							#pragma omp atomic
							next_dispersal_lattice[i1][j2][this_species - 1] += this_fecundity * distance_probability[k][l - 1] / distance_probability_sum;
							// lower left quadrant
							#pragma omp atomic
							next_dispersal_lattice[i2][j2][this_species - 1] += this_fecundity * distance_probability[k][l - 1] / distance_probability_sum;
						}
					}
					delete[] distance_probability[k];
				}
				delete[] distance_probability;
			}
			delete[] neighborhood_abundance;
		}
		else {
			// focal individual dies
			next_lattice[i][j] = 0;
			#pragma omp atomic
			species_abundance[this_species - 1]--;
		}
	}

	// no matter what happened in this site, a total of four random numbers will be discarded (the maximum number of random numbers used in the simulation)
	discardRandom(((unsigned long long) 4 - (this_random_count - start_random_count)), this_random_count, this_random_generator);
	return;
}


