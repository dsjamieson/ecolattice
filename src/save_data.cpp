
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2020			
	 *
	 *	 methods for Ecolattice class. save competition file after
	 *	 initialization. save lattice and dispersal files after 
	 *	 each time step.
	 *
	 ***********************************************************/

#include "ecolattice.hpp"

void Ecolattice::saveCompetition(void) {
	/* before the first time step, saves the parameters from this simulation to file, including species
	occupancy, survival, maximum competition, dispersal probability, dispersal length, intrinsic fecundity
	fecundity and growth competition matrices. also saves metrics calculated from the competition matrices,
	i.e., transitivity. */

	int i, j;

	std::ofstream competition_file;
	competition_file.open(outfile_base + "_competition.csv", std::ios::out | std::ios::trunc);

	if (!competition_file.is_open()) {
			if (id == 0)
				fprintf(stderr, "Error, could not open competition file to save\n");
			exit(0);
	}

	competition_file << "# Seeds:" << std::endl;
	for (i = 0; i < 5; i++) {
		competition_file << " " << seeds[i];
		if (i < 4)
			competition_file << ",";
	} 
	competition_file << std::endl;

	competition_file << "# Species Occupancy:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << species_occupancy[i];
		if (i < num_species-1)
			competition_file << ",";
	}
	competition_file << std::endl;

	competition_file << "# Juvenile Survival:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << juvenile_survival_probability[i];
		if (i < num_species-1)
			competition_file << ",";

	}
	competition_file << std::endl;

	competition_file << "# Adult Survival:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << adult_survival_probability[i];
		if (i < num_species-1)
			competition_file << ",";

	}
	competition_file << std::endl;

	competition_file << "# Maximum Competition:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << maximum_competition[i];
		if (i < num_species-1)
			competition_file << ",";
	}
	competition_file << std::endl;

	competition_file << "# Dispersal probability:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << dispersal_probability[i];
		if (i < num_species-1)
			competition_file << ",";
	}
	competition_file << std::endl;

	competition_file << "# Dispersal length:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << dispersal_length[i];
		if (i < num_species-1)
			competition_file << ",";
	}
	competition_file << std::endl;

	competition_file << "# Intrinsic fecundity:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << intrinsic_fecundity[i];
		if (i < num_species-1)
			competition_file << ",";
	}
	competition_file << std::endl;
	
	competition_file << "# Fecundity competition:" << std::endl;
	for (i = 0; i < num_species; i++) {
		for (j = 0; j < num_species; j++) {
			if (competition_fecundity[i][j] < 0)
				competition_file << competition_fecundity[i][j];
			else
				competition_file << " " << competition_fecundity[i][j];
			if (j != num_species - 1)
				competition_file << ", ";
			}
			competition_file << std::endl;
		}

	competition_file << "# Growth competition:" << std::endl;
	for (i = 0; i < num_species; i++) {
		for (j = 0; j < num_species; j++) {
			if(competition_growth[i][j] < 0)
				competition_file << competition_growth[i][j];
			else
				competition_file << " " << competition_growth[i][j];
			if (j != num_species - 1)
				competition_file << ", ";
			}
			competition_file << std::endl;
		}

	competition_file << "# Fecundity transitivity ranks:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << fecundity_rank[i];
		if (i < num_species-1)
			competition_file << ",";
	}
	competition_file << std::endl;

	competition_file << "# Growth transitivity ranks:" << std::endl;
	for (i = 0; i < num_species; i++) {
		competition_file << " " << growth_rank[i];
		if (i < num_species-1)
			competition_file << ",";
	}
	competition_file << std::endl;

	competition_file << "# Fecundity transitivity:" << std::endl;
	for (i = 0; i < num_species; i++) {
		for (j = 0; j < num_species; j++) {
			if (fecundity_transitivity[i][j] < 0)
				competition_file << fecundity_transitivity[i][j];
			else
				competition_file << " " << fecundity_transitivity[i][j];
			if (j != num_species - 1)
				competition_file << ", ";
			}
			competition_file << std::endl;
		}

	competition_file << "# Growth transitivity:" << std::endl;
	for (i = 0; i < num_species; i++) {
		for (j = 0; j < num_species; j++) {
			if (growth_transitivity[i][j] < 0)
				competition_file << growth_transitivity[i][j];
			else
				competition_file << " " << growth_transitivity[i][j];
			if (j != num_species - 1)
				competition_file << ", ";
			}
			competition_file << std::endl;
		}

	competition_file <<  "# Fecundity imbalance mean:" << std::endl  << fecundity_imbalance_mean <<  std::endl;
	competition_file <<  "# Growth imbalance mean:" << std::endl  << growth_imbalance_mean <<  std::endl;
	competition_file <<  "# Fecundity relative intransitivity:" << std::endl  << fecundity_relative_intransitivity <<  std::endl;
	competition_file <<  "# Growth relative intransitivity:" << std::endl  << growth_relative_intransitivity <<  std::endl;
	competition_file <<  "# Fecundity-growth cross correlation:" << std::endl  << fecundity_growth_correlation <<  std::endl;
	competition_file <<  "# Number of restarts:" << std::endl  << num_restarts <<  std::endl;
	
	competition_file.close();

	return;
}

void Ecolattice::saveLattice(int time_step) {
	/* saves the species locations in the lattice from the current time step to file */

	std::ofstream lattice_file;
	lattice_file.open(outfile_base+"_" + std::to_string(time_step) + ".csv", std::ios::out | std::ios::trunc);

	if (!lattice_file.is_open()) {
			if (id == 0)
				fprintf(stderr, "Error, could not open time step %d lattice file file to load\n", time_step);
			exit(0);
	}
	for (int i = 0; i < lattice_size; i++) {
		for (int j = 0; j < lattice_size; j++) {
			if (lattice[i][j] < 0)
				lattice_file << lattice[i][j];
			else
				lattice_file << " " << lattice[i][j];
			if (j != lattice_size - 1)
				lattice_file << ", ";
			}
			lattice_file << std::endl;
		}

	lattice_file.close();

	return;
}

void Ecolattice::saveDispersal(int time_step) {
	/* the dispersal lattice (locations of seeds of each species) is saved for the current time step.
	even and odd time steps are saved, so there is always one complete dispersal table to use to
	restart the simulation. */

	#pragma omp for nowait
	for (int k = 1; k < num_species + 1; k++) {
		FILE * dispersal_file;
		std::string filename;	
		filename = outfile_base + "_dispersal_s" + std::to_string(k) + "_" + std::to_string(time_step % 2) + ".csv";	
		dispersal_file = fopen(filename.c_str(), "w+");
		if (!dispersal_file) {
			if (id == 0)
				fprintf(stderr, "Error, could not open dispersal file for species %d to save\n", k);
			exit(0);
		}
		fprintf(dispersal_file, "# Dispersal for species %d, time step %d\n", k, time_step);
		for (int i = 0; i < lattice_size; i++) {
			for (int j = 0; j < lattice_size; j++) {
				fprintf(dispersal_file, "\t%.5f", dispersal_lattice[i][j][k - 1]);
				if (j < lattice_size - 1)
					fprintf(dispersal_file, ",");
			}
			fprintf(dispersal_file, "\n");
		}
		fclose(dispersal_file);
	}
	return;
}

