
	 /**********************************************************
	 * ecolattice
	 *			D.S Jamieson and N.L Kinlock, 2018			
	 *
	 *	 this file include methods for class Simulation.
	 *	 these methods save the lattice (species locations)
	 *	 and dispersal lattice (seed locations) for the current
	 *	 time step. this is so that a failed simulation
	 *	 does not have to start from the beginning.
	 *
	 ***********************************************************/

#include "simulation.h"

void Simulation::saveDispersal(int time_step) {
	/* the dispersal lattice (locations of seeds of each species) is saved for the current time step */

	int i, j, k;
	std::ofstream dispersal_file;

	for (k = 1; k < num_species + 1; k++) {

		dispersal_file.open(outfile_base + "_dispersal_s" + std::to_string(k) + "_" + std::to_string(time_step % 2) + ".csv", std::ios::out | std::ios::trunc);

		if (!dispersal_file.is_open()) {
			if (id == 0)
				fprintf(stderr, "Error, could not open dispersal file for species %d to save\n", k);
			exit(0);
		}

		dispersal_file << "# Dispersal for species " << k + 1 << ", time step " << time_step <<  std::endl;

		for (i = 0; i < lattice_size; i++) {
			for (j = 0; j < lattice_size; j++) {
				dispersal_file << dispersal_lattice[i][j][k - 1] << ", ";
			}
			dispersal_file << std::endl;
		}
		dispersal_file.close();
	}

	return;

}


void Simulation::loadLattice() {
	/* load the lattice, or the species locations, from the previous simulation. */

	int i, j, k, col_num;
	int line_num = 0;
	std::ifstream lattice_file;
	std::string line;
	std::string value;

	lattice_file.open(outfile_base + "_" + std::to_string(restart_time) + ".csv", std::ios::out | std::ios::trunc);

	if (!lattice_file.is_open()) {
		if (id == 0)
			fprintf(stderr, "Error, could not open lattice file for time step %d to load\n", restart_time);
		exit(0);
	}

	while (getline(lattice_file, line)) {
		line_num++;
		if (line_num > lattice_size) {
			if (id == 0)
				fprintf(stderr, "Error, too many lines in lattice file for time step %d\n", restart_time);
			exit(0);
		}
		
		col_num = 0;

		std::istringstream values(line);

		while (values >> value) {
			col_num++;
			if (col_num > lattice_size) {
				if (id == 0)
					fprintf(stderr, "Error, too many columns in lattice file for time step %d, line %d\n", restart_time, line_num);
				exit(0);
			}

			value = value.substr(0, value.size() - 1);

			if (value.find_first_not_of("0123456789") == std::string::npos) {
					try {
						lattice[line_num - 1][col_num - 1] = stoi(value);
					}
					catch (...) {
						if (id == 0)
							fprintf(stderr, "Error, could not convert lattice file value given for time step %d, in line %d column %d, to positive integers\n", restart_time, line_num, col_num);
						exit(0);
					}
			}
			else {
				if (id == 0)
					fprintf(stderr, "Error, invalid lattice file value given for time step %d, in line %d column %d, to positive integers\n", restart_time, line_num, col_num );
				exit(0);
			}
		}
		if (col_num < lattice_size) {
			if (id == 0)
				fprintf(stderr, "Error, not enough columns in lattice file for time step %d, line %d\n", restart_time, line_num );
			exit(0);
		}
	}
		if (line_num > lattice_size) {
			if (id == 0)
				fprintf(stderr, "Error, not enough lines in lattice file for time step %d\n", restart_time);
			exit(0);
		}

	lattice_file.close();

	return;
}



void Simulation::loadDispersal() {
	/* load the dispersal lattice, or the seed locations, from the previous simulation. */

	int i, j, k, col_num;
	int line_num = 0;
	std::ifstream dispersal_file;
	std::string line;
	std::string value;

	for (k = 1; k < num_species + 1; k++) {

		dispersal_file.open(outfile_base + "_dispersal_s" + std::to_string(k) + "_" + std::to_string(restart_time % 2) + ".csv", std::ios::out | std::ios::trunc);

		if (!dispersal_file.is_open()) {
			if (id == 0)
				fprintf(stderr, "Error, could not open dispersal file for species %d, to load\n", k);
			exit(0);
		}
	
		getline(dispersal_file, line);
		if (trimString(line).size() != 0 ) {
			if (id == 0)
				fprintf(stderr, "Error, dispersal file for species %d doesn't begin with comment line\n", k);
			exit(0);
		}
		while (getline(dispersal_file, line)) {
			line_num++;
			if (line_num > lattice_size) {
				if (id == 0)
					fprintf(stderr, "Error, too many lines in dispersal file for species %d\n", k);
				exit(0);
			}
			
			col_num = 0;

			std::istringstream values(line);

			while (values >> value) {

				col_num++;
				if (col_num > lattice_size) {
					if (id == 0)
						fprintf(stderr, "Error, too many columns in dispersal file for species %d, line %d\n", k, line_num );
					exit(0);
				}

				value = value.substr(0, value.size() - 1);

				if (value.find_first_not_of("0123456789.") == std::string::npos) {
						try {
							dispersal_lattice[line_num - 1][col_num - 1][k - 1] = stod(value);
						}
						catch (...) {
							if (id == 0)
								fprintf(stderr, "Error, could not convert value given for species %d, in line %d column %d, to double\n", k + 1, line_num, col_num);
							exit(0);
						}
				}
				else {
					if (id == 0)
						fprintf(stderr, "Error, invalid dispersal file value given for species %d, in line %d column %d, to double\n", k + 1, line_num, col_num );
					exit(0);
				}
			}
			if (col_num < lattice_size) {
				if (id == 0)
					fprintf(stderr, "Error, not enough columns in dispersal file for species %d, line %d\n", k, line_num );
				exit(0);
			}
		}
			if (line_num > lattice_size) {
				if (id == 0)
					fprintf(stderr, "Error, not enough lines in dispersal file for species %d\n", k);
				exit(0);
			}

		dispersal_file.close();

	}

	return;

}


void Simulation::loadCompetition() {
	/* load the growth and fecundity competition matrices as well as other important parameters, from the previous simulation. */

	int i, j, col_num;
	int line_num = 0;
	std::string line;
	std::string value;
	std::string name;

	const char *names[] = {  "SpeciesOccupancy", "JuvenileSurvival",  "AdultSurvival",  "MaximumCompetition",
						     "DispersalProbability", "DispersalLength", "Fecundity" };

	double **parameters = new double*[7];

	parameters[0] = species_occupancy;
	parameters[1] = juvenile_survival_probability;
	parameters[2] = adult_survival_probability;
	parameters[3] = maximum_competition;
	parameters[4] = dispersal_probability;
	parameters[5] = dispersal_length;
	parameters[6] = intrinsic_fecundity;

	std::ifstream competition_file;
	competition_file.open(competition_filename);

	if (!competition_file.is_open()) {
			if (id == 0)
				fprintf(stderr, "Error, could not open competition file to load\n");
			exit(0);
	}

	while (getline(competition_file, line)) {
		line_num++;
		col_num = 0;
		if (trimString(line).size() != 0) {
				if (id == 0)
					fprintf(stderr, "Error, formatting issue, line %d in competition file is a not comment\n", line_num);
				exit(0);
		}

		getline(competition_file, line);
		line_num++;

		std::istringstream values( line );

		while (values >> value) {

			col_num++;
			if (col_num > num_species) {
				if (id == 0)
					fprintf(stderr, "Error, too many columns in competition file, line %d\n", line_num );
				exit(0);
			}

			if (col_num == 1) {
				if (value.compare("reset") == 0 && line_num / 2 - 1 < 5 ) {
					std::string name(names[ line_num/2 - 1]);
					int type = 2;
					if (line_num == 2)
						type = 3;
					if (line_num/2 - 1 > 4)
						type = 4;
					setRandomParameter(parameters[line_num / 2 - 1], num_species, name, type);
					break;
				}
			}

			if(col_num < num_species)
				value = value.substr(0, value.size() - 1);

			if (value.find_first_not_of("0123456789.") == std::string::npos) {

				try {
					parameters[line_num / 2 - 1][col_num - 1] = stod(value);
					}
				catch (...) {
					if (id == 0)
						fprintf(stderr, "Error, could not convert value given for competition file line %d column %d to double\n", line_num, col_num );
					exit(0);
				}
			}
			else {
				if (id == 0)
					fprintf(stderr, "Error, invalid competition file value given in line %d column %d, to double\n", line_num, col_num);
				exit(0);
			}

		}

		if (line_num == 14)
			break;

	}

	if (line_num < 14) {
		if (id == 0)
			fprintf(stderr, "Error, not enough lines from species specific parameters in competition file\n");
		exit(0);

	}

	if (getline(competition_file, line)) {
		line_num++;
		if (trimString(line).size() != 0 ) {
				if (id == 0)
					fprintf(stderr, "Error, formatting issue, line %d in competition file is a not comment\n", line_num);
				exit(0);
		}
	}
	else {
		if (id == 0)
			fprintf(stderr, "Error, competition file ends before fecundity competition loaded\n");
		exit(0);
	}

	while (getline(competition_file, line)) {

		line_num++;
		col_num = 0;

		std::istringstream values( line );

			while (values >> value) {

				col_num++;
				if (col_num > num_species) {
					if (id == 0)
						fprintf(stderr, "Error, too many columns in competition file, line %d\n", line_num);
					exit(0);
				}			

				if(col_num < num_species)
					value = value.substr(0, value.size() - 1);

				if (value.find_first_not_of("-0123456789.") == std::string::npos) {
					try {
						competition_fecundity[line_num - 16][col_num - 1] = stod(value);
					}
					catch (...) {
						if (id == 0)
							fprintf(stderr, "Error, could not convert value given for competition file line %d column %d to double\n", line_num, col_num);
						exit(0);
					}
				}
				else {
					if (id == 0)
						fprintf(stderr, "Error, invalid competition file value given in line %d column %d, to double\n", line_num, col_num);
					exit(0);
				}
			}
		if (line_num == 15 + num_species)
			break;
	}
	if (line_num < 15 + num_species) {
		if (id == 0)
			fprintf(stderr, "Error, not enough lines from fecundity competition in competition file, given the number of species\n");
		exit(0);
	}
	if (getline(competition_file, line)) {
		line_num++;
		if (trimString(line).size() != 0) {
				if (id == 0)
					fprintf(stderr, "Error, formatting issue, line %d in competition file is a not comment\n", line_num);
				exit(0);
		}
	}
	else {
		if (id == 0)
			fprintf(stderr, "Error, competition file ends before growth competition loaded\n");
		exit(0);
	}

	while (getline(competition_file, line)) {
		line_num++;
		col_num = 0;

		std::istringstream values( line );

			while (values >> value) {
				col_num++;
				if (col_num > num_species) {
					if (id == 0)
						fprintf(stderr, "Error, too many columns in competition file, line %d\n", line_num);
					exit(0);
				}

				if(col_num < num_species)
					value = value.substr(0, value.size() - 1);

				if (value.find_first_not_of("-0123456789.") == std::string::npos) {

						try {
							competition_growth[line_num - 17 - num_species][col_num - 1] = stod(value);
						}
						catch (...) {
							if (id == 0)
								fprintf(stderr, "Error, could not convert value given for competition file line %d column %d to double\n", line_num, col_num);
							exit(0);
						}
				}
				else {
					if (id == 0)
						fprintf(stderr, "Error, invalid competition file value given in line %d column %d, to double\n", line_num, col_num);
					exit(0);
				}
			}
		if (line_num == 16 + 2 * num_species)
			break;
	}
	if(line_num < 16 + 2 * num_species) {
		if (id == 0)
			fprintf(stderr, "Error, not enough lines from growth competition in competition file, given the number of species\n");
		exit(0);

	}

	delete[] parameters;

	competition_file.close();

	return;

}



int Simulation::getRestartTime() {
	return restart_time;
}
