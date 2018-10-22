
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

	int numprocs, myid;

   	MPI_Init(&argc,&argv);                 // Initialize
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);  // Get # processors
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);      // Get my rank (id)
	MPI_Status status; 

	// Multiproccessor MPI version
	if(numprocs > 1) {

		Simulation sim(argv[1], myid);

		int i,j,k,l,m, box_size, total_workers, num_species, time_step, max_time_step;
		int this_site, buffer_size, this_species;
		int *task_site_number;
		double *buffer;

		box_size = sim.getBoxSize();
		num_species = sim.getSpecies();
		max_time_step = sim.getMaxTimeStep();

		// Allocate buffer
		buffer_size = box_size*box_size*(1+num_species); 
		buffer = new double[buffer_size];
		if(!buffer) {
			fprintf(stderr, "Error, could not allocate memory for buffer in main\n");
			MPI_Finalize();
			exit(-1);
		}
		for(i=0; i<buffer_size; i++)
			buffer[i] = 0.;

		// Assign sites to processors
		task_site_number = new int[numprocs];
		if(!task_site_number) {
			fprintf(stderr, "Error, could not allocate memory for parallel tasks in main\n");
			MPI_Finalize();
			exit(-1);
		}
		task_site_number[0] = 0;
		for(i=1;i<numprocs;i++) {

			task_site_number[i] = task_site_number[i-1] + box_size*box_size/(numprocs-1);
			if( i <= box_size*box_size%(numprocs-1) )
				task_site_number[i]++;
			if( task_site_number[i] <= box_size*box_size  )
				total_workers = i;
		}

		// Central task
		if(myid == 0) {

			// Save initial data
			sim.saveCompetition();
			sim.saveProperties();
			sim.saveBox(0);
			sim.resetBox();

			for(time_step=1;time_step<max_time_step+1;time_step++) {

				// Listen to workers for time step output
				for(m=0; m<total_workers; m++) {

					MPI_Recv(buffer, buffer_size, MPI_DOUBLE, MPI_ANY_SOURCE, time_step, MPI_COMM_WORLD, &status);

					// Update dispersal for this species
					for(i=0; i<box_size; i++) {
						for(j=0; j<box_size; j++ ) {

							// Set box element i,j
							sim.addSite(i,j, (int) buffer[j + i*box_size] );
		
							for(k=0; k<num_species; k++) 
								sim.addDispersal(i,j,k, buffer[ box_size*box_size + j + i*box_size + k*box_size*box_size ]);		
						}
					}
				}

				// Fill buffer for sending			
				for(i=0; i<box_size;i++) {
					for(j=0; j<box_size; j++) {

						buffer[j+i*box_size] = (double) sim.getSite(i,j);

						for(k=0;k<num_species;k++) {
							buffer[ box_size*box_size + j + i*box_size + k*box_size*box_size  ] = sim.getDispersal(i,j,k);
						}
					}
				}

				// Send updated box and disperal to workers for next timestep
				for(k=0; k<total_workers; k++)
					MPI_Send(buffer, buffer_size, MPI_DOUBLE, k+1, time_step, MPI_COMM_WORLD);

				// Save updated box and reset box and dispersal
				sim.saveBox(time_step);
				sim.resetBox();

				fprintf(stdout, "done step: %d\n", time_step);

			}

		}
		// Workers
		else if( task_site_number[myid-1] < box_size*box_size  ) {

			for( time_step=1; time_step<max_time_step+1; time_step++) {

				// Clear buffer
				for(k=0;k<buffer_size;k++)
					buffer[k] = 0.;

				for(this_site=task_site_number[myid-1]; this_site<task_site_number[myid]; this_site++) {

					i = this_site/box_size;
					j = this_site - i*box_size;

					sim.discardRandom( (  ( (unsigned long long) ( 4*box_size*box_size*(time_step-1) + 4*(j+box_size*i) ) ) - ( sim.getRandomCount() ) ) );

					sim.updateSingleSite(i,j);

				}


				for(this_site=task_site_number[myid-1]; this_site<task_site_number[myid]; this_site++) {

					i = this_site/box_size;
					j = this_site - i*box_size;

					// Fill buffer for sending to central process
					buffer[j+box_size*i] = (double) sim.getNextSite(i,j);
					this_species = abs( sim.getSite(i,j) );

					if(this_species!=0) {
						for(k=0; k<box_size; k++) {
							for(l=0; l<box_size; l++) {

								buffer[box_size*box_size + l + k*box_size + (this_species-1)*box_size*box_size ] = sim.getNextDispersal(k, l, this_species-1 );

							}
						}
					}
				}

				// Send i, j timestep results 
				MPI_Send(buffer, buffer_size, MPI_DOUBLE, 0, time_step, MPI_COMM_WORLD);

				// Receive update box and dispersal
				MPI_Recv(buffer, buffer_size, MPI_DOUBLE, 0, time_step, MPI_COMM_WORLD, &status);

				// Update box and dispersal for next time step
				if(time_step!=max_time_step) {
					for(i=0; i<box_size; i++) {
						for(j=0; j<box_size; j++) {

							sim.setSite(i,j, (int) buffer[ j + i*box_size ]);

							for(k=0; k<num_species; k++) {
								sim.setDispersal(i,j,k, buffer[box_size*box_size + j + i*box_size + k*box_size*box_size]);
							}
						}
					}			
				}

				sim.resetNextBox();

			}
		}

		delete[] buffer;
		delete[] task_site_number;

		MPI_Finalize();
		return(0);

	}
	// Single processor version
	else {

		int i,j, box_size, time_step, max_time_step;

		Simulation sim(argv[1], myid);
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

		MPI_Finalize();
		return(0);

	}
}


