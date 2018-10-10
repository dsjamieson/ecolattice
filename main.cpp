#include <iostream>
#include <random>
#include <vector>
//#include <omp.h>
#include "simulation.h"

using namespace std;
mt19937 global_random_generator;

int main(int argc, char* argv[]) {

	if(argc != 2){
		fprintf(stderr, "Usage: ecolat input_filename\n");
		exit(0);
	}

	int i,j, box_size, time_step, max_time_step;

	Simulation sim(argv[1]);
	sim.saveBox(0);

	box_size = sim.getBoxSize();
	max_time_step = sim.getMaxTimeStep();
	int new_seed = sim.getNewSeed(1);
	mt19937 local_random_generator(new_seed);	

	for(time_step=1;time_step<max_time_step+1;time_step++) {

		for(i=0; i<box_size; i++) {
			for(j=0;j<box_size;j++) {
				sim.updateSingleSite(i,j,local_random_generator);
			}
		}
	
		sim.nextToThis();
		sim.saveBox(time_step);

	}

	return(0);
}


