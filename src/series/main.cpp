
	 /**********************************************************
	 *
	 *			D.S Jamieson and N.L Kinlock, 2018			
	 *
	 *	 
	 *	 
	 *	
	 *	 
	 *	
	 *
	 ***********************************************************/

#include "simulation.h"

int main(int argc, char* argv[]) {

	if(argc != 2){
		fprintf(stderr, "Usage: ecolat input_filename\n");
		exit(0);
	}

	int i,j, box_size, time_step, max_time_step;

	Simulation sim(argv[1], 0);
	sim.saveCompetition();
	sim.saveProperties();
	sim.saveBox(0);

	box_size = sim.getBoxSize();
	max_time_step = sim.getMaxTimeStep();

	for(time_step=1;time_step<max_time_step+1;time_step++) {

		for(i=0; i<box_size; i++) {
			for(j=0;j<box_size;j++) {

				sim.discardRandom( ( (int64_t)  ( 4*box_size*box_size*(time_step-1) + 4*(j+box_size*i) ) ) - sim.getRandomCount() );
				sim.updateSingleSite(i,j);

			}
		}
	
		sim.nextToThis();
		sim.saveBox(time_step);

	}

	return(0);
}


