
	 /**********************************************************
	 * ecolattice
	 *			D.S Jamieson and N.L Kinlock, 2020			
	 *
	 *	 methods for Ecolattice class. random draws using Mersenne
	 *	 Twister RNG, aligned so simulations are repeatable.
	 *
	 ***********************************************************/

#include "ecolattice.hpp"

void Ecolattice::drawRandomSeeds(void) {
	/* draw random seeds from random device (5) and system clock (1). system clock used for systems 
	that do not have random device capability. */
	int i;
	std::random_device r;
	seeds[0] = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	for (i = 1; i < 5; i++)
		seeds[i] = (unsigned int) r(); 

	return;
}

void Ecolattice::seedRandomGenerator(void) {
	/* create a random vector of seeds (a seed sequence) given the seeds specified randomly or in
	the parameter file (restart simulation). seeds fed to the global RNG. */
    std::seed_seq seq(seeds, seeds + 5);
	std::vector<std::uint32_t> seed_vector(std::mt19937::state_size);
    seq.generate(seed_vector.begin(), seed_vector.end());
	std::seed_seq seq2(seed_vector.begin(), seed_vector.end());
	global_random_generator.seed(seq2);
	return;
}

void Ecolattice::discardRandom(unsigned long long t_num_discard) {
	/* add to the random count and discard values from the global RNG. */
	global_random_generator.discard(t_num_discard);
	random_count += t_num_discard;
	return;
}

std::mt19937& Ecolattice::generateRandom(void) {
	/* add to the random count and get a random draw from the global RNG. */
	random_count += 2;
	return global_random_generator;
}

unsigned int Ecolattice::getRandom(void) {
	random_count++;
	return global_random_generator();

}

double Ecolattice::getRandomUniformReal(double t_lower_bound, double t_upper_bound) {
	return (double) t_lower_bound + (t_upper_bound - t_lower_bound) * getRandom() / static_cast<double>(std::mt19937::max()) ;
}

double Ecolattice::getRandomNormal(double t_mean, double t_sdev) {
	static bool hasSpare = false;
	static double spare;
	double s = 0, u = 0, v = 0; 
	if(hasSpare) {
		hasSpare = false;
		return t_mean + t_sdev * spare;
	}
	hasSpare = true;
	while( (s >= 1.0) || (s == 0) ) {
		u = (double) getRandom() / static_cast<double>(std::mt19937::max()) * 2.0 - 1.0;
		v = (double) getRandom() / static_cast<double>(std::mt19937::max()) * 2.0 - 1.0;
		s = u * u + v * v;
	}
	s = sqrt(-2.0 * log(s) / s);
	spare = v * s;
	return t_mean + t_sdev * u * s;
}

void Ecolattice::initializeRandomParameter(std::vector<double> & parameter_value, std::string parameter_name, int type) {
	/* some parameters in the data file  are specified as parameters of a normal distribution (e.g., intrinsic fecundity). 
	this method initializes the mean and standard deviation arrays to be sent to the method 'initializeNormalArray' which 
	will initialize the actual array of random draws. this method is specific to parameters that are probabilities. */

	// type = 0 -> not essential, does not need to be set in parameter file
	// type = 1 -> essential, must be set in parameter file
	// type = 2 -> essential, all values must be between 0 and 1
	// type = 3 -> essential, treated as weights
	// type = 4 -> essential, must be non-negative

	std::vector<double> mean(num_species, 0.);
	std::vector<double> sdev(num_species, 0.);
	getParameter(mean, parameter_name, type);
	getParameter(sdev, parameter_name + "Sdev", 0);

	initializeNormalRandomArray(parameter_value, mean, sdev);

	return;

}

void Ecolattice::initializeNormalRandomArray(std::vector<double> & array, std::vector<double> & mean, std::vector<double> & sdev) {
	/* for parameters that are parameters of a normal distribution (e.g. intrinsic fecundity). 
	this method takes draws from normal distributions defined by the mean and standard deviation arrays. */

	unsigned long i;
	for (i = 0; i < array.size(); i++)
		array[i] = fabs(getRandomNormal(mean[i], sdev[i]));

	return;

}
