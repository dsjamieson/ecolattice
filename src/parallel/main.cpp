
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2018			
	 *
	 *	 
	 *	 this is the main body of the simulation program.
	 *	 the program can run in parallel, using MPI, or
	 *	 in serial, depending on the number of processors.
	 *	 MPI must be installed on your computer to run this
	 *	 version.
	 *
	 ***********************************************************/

#include "simulation.h"

int main(int argc, char* argv[]) {

	if (argc != 2) {
		fprintf(stderr, "Usage: ecolattice input_filename\n");
		exit(0);
	}

	int numprocs, myid;

   	MPI_Init(&argc, &argv);   // initialize MPI
    	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);  // get # processors
    	MPI_Comm_rank(MPI_COMM_WORLD, &myid);  // get rank (id)
	MPI_Status status; 

	// multiproccessor version
	if (numprocs > 1) {
		
		// instantiating class Simulation with object sim
		// requires file name with list of parameters and worker ID number
		Simulation sim(argv[1], myid);

		int i, j, k, l, m, lattice_size, total_workers, num_species, persistence;
		int start_time, time_step, max_time_step, min_persistence;
		int this_site, buffer_size, this_species, this_stage;
		int *task_site_number, *species_extant;
		double *buffer;

		lattice_size = sim.getLatticeSize();  // number of cells in one dimension of the lattice
		num_species = sim.getSpecies();  // number of species in the simulation
		start_time = sim.getRestartTime() + 1;
		min_persistence = sim.getMinPersistence();
		max_time_step = sim.getMaxTimeStep();  // duration of the simulation
		for ( i = 0; i < num_species; i++)
			species_extant[i] = 0;
		persistence = 0;


		// allocate buffer
		// buffer is an array that stores the lattice grid, with species locations, and the dispersal lattice grid, with seed locations
		// buffer is used to pass information between central task and workers
		buffer_size = lattice_size * lattice_size * (1 + num_species) + num_species; 
		buffer = new double[buffer_size];
		if (!buffer) {
			fprintf(stderr, "Error, could not allocate memory for buffer in 'main' function\n");
			MPI_Finalize();
			exit(-1);
		}
		for (i = 0; i < buffer_size; i++)
			buffer[i] = 0.;

		// assign sites to processors
		// lattice is split among the processors to parallelize each step of the simulation.
		// each processor starts at evenly spaced sections of the lattice, split by cells read left to right (not in subgrids)
		task_site_number = new int[numprocs];
		if (!task_site_number) {
			fprintf(stderr, "Error, could not allocate memory for parallel tasks in 'main' function\n");
			MPI_Finalize();
			exit(-1);
		}
		task_site_number[0] = 0;
		for (i = 1; i < numprocs; i++) {
			task_site_number[i] = task_site_number[i - 1] + lattice_size * lattice_size / (numprocs - 1);
			if (i <= lattice_size * lattice_size % (numprocs - 1))
				task_site_number[i]++;
			if (task_site_number[i] <= lattice_size * lattice_size)
				total_workers = i;
		}

		// central task
		// central task does not have lattice sites assigned to it. it receives the partially complete buffer from each worker and sends the complete buffer to each worker.
		if (myid == 0) {

			// central task saves initial state at time 0, including parameters, competition matrices, and initial individual locations
			if (start_time == 1) {
				sim.saveCompetition();
				sim.saveLattice(0);
				sim.resetLattice();
			}

			for (time_step = start_time; time_step < max_time_step + 1; time_step++) {

				// central task receives partially complete buffer from each worker, including species and seed positions at sites assigned to that worker
				// adds these locations to the array 'lattice' as a part of the object 'sim' (and the same for dispersal lattice)
				for (m = 0; m < total_workers; m++) {

					MPI_Recv(buffer, buffer_size, MPI_DOUBLE, MPI_ANY_SOURCE, time_step, MPI_COMM_WORLD, &status);

					for (i = 0; i < lattice_size; i++) {
						for (j = 0; j < lattice_size; j++) {
							sim.addSite(i, j, (int) buffer[j + i * lattice_size]);
							for (k = 0; k < num_species; k++) 
								sim.addDispersal(i, j, k, buffer[lattice_size * lattice_size + j + i * lattice_size + k * lattice_size * lattice_size]);		
						}
					}

					for (i = 0; i < num_species; i++) {
						if (buffer[buffer_size - num_species + i] != 0)
							species_extant[i] = 1;
					}

				}

				for (i = 0; i < num_species; i++) {
					if (species_extant[i] != 0)
						persistence++;
				}

				if (persistence < min_persistence) {
			
					sim.setRandomSeeds();
					for (i = 0; i < 5; i++)
						buffer[i] = (double) sim.getSeed(i); 	

					for (k = 0; k < total_workers; k++)
						MPI_Send(buffer, buffer_size, MPI_DOUBLE, k+1, -1, MPI_COMM_WORLD);
					sim.reinitializeSim(time_step);
					sim.saveCompetition();
					sim.saveLattice(0);
					sim.resetLattice();
					for ( i = 0; i < num_species; i++)
						species_extant[i] = 0;
					persistence = 0;
					time_step = 0;					
					continue;

				}

				// central task updates its own buffer to completeness, including all species and seed locations	
				for (i = 0; i < lattice_size; i++) {
					for (j = 0; j < lattice_size; j++) {
						buffer[j + i * lattice_size] = (double) sim.getSite(i, j);
						for (k = 0; k < num_species; k++) {
							buffer[lattice_size * lattice_size + j + i * lattice_size + k * lattice_size * lattice_size] = sim.getDispersal(i, j, k);
						}
					}
				}

				// central task sends complete buffer to workers
				for (k = 0; k < total_workers; k++)
					MPI_Send(buffer, buffer_size, MPI_DOUBLE, k+1, time_step, MPI_COMM_WORLD);

				// central task saves lattice from current time step to file and resets lattice and dispersal lattice for next time step
				sim.saveLattice(time_step);
				sim.saveDispersal(time_step);
				sim.resetLattice();
				for ( i = 0; i < num_species; i++)
					species_extant[i] = 0;
				persistence = 0;

				fprintf(stdout, "done step: %d\n", time_step);
			}
		}

		// workers
		// each worker is assigned a certain number of sites at which it will update the species and seed locations in 'next_lattice' and 'next_dispersal_lattice'
		else if (task_site_number[myid - 1] < lattice_size * lattice_size) {

			for (time_step = start_time; time_step < max_time_step + 1; time_step++) {

				// worker clears its copy of buffer before updating sites
				for (k = 0; k < buffer_size; k++)
					buffer[k] = 0.;

				// worker is assigned sites i, j given ID number
				for (this_site = task_site_number[myid - 1]; this_site < task_site_number[myid]; this_site++) {
					i = this_site / lattice_size;
					j = this_site - i * lattice_size;

					// RNG discards values so that the simulations are repeatable
					sim.discardRandom((((unsigned long long) (4 * lattice_size * lattice_size * (time_step - 1) + 4 * (j + lattice_size * i))) - (sim.getRandomCount())));
					
					// worker updates single site in 'next_lattice' and 'next_dispersal_lattice'
					sim.updateSingleSite(i, j);
				}


				// worker updates its copy of buffer (partially complete) to send to central task
				for (this_site = task_site_number[myid - 1]; this_site < task_site_number[myid]; this_site++) {

					i = this_site / lattice_size;
					j = this_site - i * lattice_size;

					// worker updates its copy of buffer with species and seed locations from 'next_lattice' and 'next_dispersal_lattice'

					m = sim.getNextSite(i, j)
					buffer[j + lattice_size * i] = (double) m;
					if( m != 0 && species_extant[m - 1] != 0)
						species_extant[m - 1] = 1;

					this_species = abs(sim.getSite(i, j));
					if (this_species != 0)
						this_stage = (int) sim.getSite(i, j) / this_species;
					else
						this_stage = 0;
					if (this_stage > 0) {
						for (k = 0; k < lattice_size; k++) {
							for (l = 0; l < lattice_size; l++) {
								buffer[lattice_size * lattice_size + l + k * lattice_size + (this_species - 1) * lattice_size * lattice_size] = sim.getNextDispersal(k, l, this_species - 1);
							}
						}
					}
				}

				for (i = 0; i < num_species; i++)
					buffer[buffer_size - num_species + i] =  species_extant[i];

				// worker sends partially complete buffer to central task
				MPI_Send(buffer, buffer_size, MPI_DOUBLE, 0, time_step, MPI_COMM_WORLD);

				// worker receives complete buffer from central task
				MPI_Recv(buffer, buffer_size, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

				if (status.MPI_Tag == -1) {

					for (i = 0; i < 5; i++) 
						sim.setSeed(i, (unsigned int) buffer[i]);

					sim.reinitializeSimulation(time_step);
					sim.resetNextLattice();
					for ( i = 0; i < num_species; i++)
						species_extant[i] = 0;
					time_step = 0;
					continue;
		
				}
				else if ( status.MPI_Tag != time_step  ) {

					fprintf(stderr, "Error, proc %d received time step %d, should be %d\n", myid, status.MPI_Tag, time_step);
					MPI_Finalize();
					exit(-1);

				}

				// worker updates 'lattice' and 'dispersal_lattice' with information from complete buffer
				// to be used in the next time step
				if (time_step != max_time_step) {
					for (i = 0; i < lattice_size; i++) {
						for (j = 0; j < lattice_size; j++) {
							sim.setSite(i, j, (int) buffer[j + i * lattice_size]);
							for (k = 0; k < num_species; k++) {
								sim.setDispersal(i, j, k, buffer[lattice_size * lattice_size + j + i * lattice_size + k * lattice_size * lattice_size]);
							}
						}
					}			
				}
				// worker clears 'next_lattice' and 'next_dispersal_lattice', which will be filled for the next time step
				sim.resetNextLattice();
				for ( i = 0; i < num_species; i++)
					species_extant[i] = 0;

			}
		}

		delete[] buffer;
		delete[] task_site_number;

	}

	// single processor version
	else {

		int i, j, lattice_size, time_step, max_time_step;

		Simulation sim(argv[1], myid);
		// save initial values, including competition matrices, initial species locations
		sim.saveCompetition();
		sim.saveProperties();
		sim.saveLattice(0);

		lattice_size = sim.getLatticeSize();  // length of one dimension of the lattice
		max_time_step = sim.getMaxTimeStep();  // duration of simulation

		for (time_step = 1; time_step < max_time_step + 1; time_step++) {
			for (i = 0; i < lattice_size; i++) {
				for (j = 0; j < lattice_size; j++) {
					// discard values from the RNG so simulation is repeatable
					sim.discardRandom(((int64_t)  (4 * lattice_size * lattice_size * (time_step - 1) + 4 * (j + lattice_size * i))) - sim.getRandomCount());
					sim.updateSingleSite(i, j);  // update current site
				}
			}
			// move to next site
			sim.nextToThis();
			// save lattice from current time step
			sim.saveLattice(time_step);
		}

	}

	MPI_Finalize();
	return(0);

}


