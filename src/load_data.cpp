
	 /**********************************************************
	 * ecolattice
	 *			D.S Jamieson and N.L Kinlock, 2018			
	 *
	 *	 methods for Ecolattice class. load the lattice 
	 *	 (species locations) and dispersal lattice (seed locations) 
	 *	 for a restart or continue simulation.
	 *
	 ***********************************************************/

#include "ecolattice.hpp"

void Ecolattice::loadSeeds(void) {
	/* if simulation in restarted, then the seeds from the previous simulation are read in from the competition file,
	in the format, "OutfileBase_competition.csv" */

	int col_num;
	std::string line;
	std::string value;

	std::ifstream competition_file;
	competition_file.open(competition_filename);

	if (!competition_file.is_open()) {
			if (id == 0)
				fprintf(stderr, "Error, could not open competition file to load\n");
			exit(0);
	}

	if (getline(competition_file, line)) {
		if (trimString(line).size() != 0 ) {
				if (id == 0)
					fprintf(stderr, "Error, formatting issue, line 1 in competition file is a not comment\n");
				exit(0);
		}
	}
	else {
		if (id == 0)
			fprintf(stderr, "Error, competition file ends before seeds are loaded\n");
		exit(0);
	}

	if (getline(competition_file, line)) {

		col_num = 0;
		std::istringstream values(line);

		while (values >> value) {
			col_num++;
			if (col_num > 5) {
				if (id == 0)
					fprintf(stderr, "Error, competition file must have 5 seeds, too many seeds given\n");
				exit(0);
			}

			if (col_num < 5)
				value = value.substr(0, value.size() - 1);

			if (value.find_first_not_of("0123456789") == std::string::npos) {

				try {
					seeds[col_num - 1] = (unsigned int) stoul(value);
					}
				catch (...) {
					if (id == 0)
						fprintf(stderr, "Error, could not convert value given for competition file column %d to int\n", col_num);
					exit(0);
				}
			}
			else {
				if (id == 0)
					fprintf(stderr, "Error, invalid competition file value given in column %d, to int\n", col_num);
				exit(0);
			}

		}

		if (col_num < 5) {
			if (id == 0)
				fprintf(stderr, "Error, competition file must have 5 seeds, too few given\n");
			exit(0);
		}
	}
	else{
		if (id == 0)
			fprintf(stderr, "Error, no seeds in competition file\n");
		exit(0);
	}
	competition_file.close();
	return;
}

void Ecolattice::loadLattice(void) {
	/* load the lattice, or the species locations, from another simulation. time step reloaded based on
	"ContinueTime." outfile base must match. */

	int col_num;
	int line_num = 0;
	std::ifstream lattice_file;
	std::string line;
	std::string value;

	lattice_file.open(outfile_base + "_" + std::to_string(continue_time) + ".csv");

	if (!lattice_file.is_open()) {
		if (id == 0)
			fprintf(stderr, "Error, could not open lattice file for time step %d to load\n", continue_time);
		exit(0);
	}

	while (getline(lattice_file, line)) {
		line_num++;
		if (line_num > lattice_size) {
			if (id == 0)
				fprintf(stderr, "Error, too many lines in lattice file for time step %d\n", continue_time);
			exit(0);
		}
		
		col_num = 0;

		std::istringstream values(line);

		while (values >> value) {
			col_num++;
			if (col_num > lattice_size) {
				if (id == 0)
					fprintf(stderr, "Error, too many columns in lattice file for time step %d, line %d\n", continue_time, line_num);
				exit(0);
			}

			if (col_num < lattice_size)
				value = value.substr(0, value.size() - 1);

			if (value.find_first_not_of("-0123456789") == std::string::npos) {
					try {
						lattice[line_num - 1][col_num - 1] = stoi(value);
					}
					catch (...) {
						if (id == 0)
							fprintf(stderr, "Error, could not convert lattice file value given for time step %d, in line %d column %d, to positive integers\n", continue_time, line_num, col_num);
						exit(0);
					}
			}
			else {
				if (id == 0)
					fprintf(stderr, "Error, invalid lattice file value given for time step %d, in line %d column %d, to positive integers\n", continue_time, line_num, col_num);
				exit(0);
			}
		}
		if (col_num < lattice_size) {
			if (id == 0)
				fprintf(stderr, "Error, not enough columns in lattice file for time step %d, line %d\n", continue_time, line_num);
			exit(0);
		}
	}
		if (line_num > lattice_size) {
			if (id == 0)
				fprintf(stderr, "Error, not enough lines in lattice file for time step %d\n", continue_time);
			exit(0);
		}

	lattice_file.close();

	return;
}

void Ecolattice::loadDispersal(void) {
	/* load the dispersal lattice, or the seed locations, from the previous simulation. time step reloaded based on
	"ContinueTime." outfile base must match. */

	int k, col_num;
	int line_num = 0;
	std::ifstream dispersal_file;
	std::string line;
	std::string value;

	for (k = 1; k < num_species + 1; k++) {

		dispersal_file.open(outfile_base + "_dispersal_s" + std::to_string(k) + "_" + std::to_string(continue_time % 2) + ".csv");
		line_num = 0;

		if (!dispersal_file.is_open()) {
			if (id == 0)
				fprintf(stderr, "Error, could not open dispersal file for species %d\n", k);
			exit(0);
		}
	
		getline(dispersal_file, line);
		if (trimString(line).size() != 0 ) {
			if (id == 0)
				fprintf(stderr, "Error, dispersal file for species %d doesn't begin with comment line\n", k);
			exit(0);
		}

		std::istringstream time_check_values(line);

		while (time_check_values >> value) 	
			continue;

		int check_continue_time = stoi(value);
		if (continue_time != check_continue_time) {
			fprintf(stderr, "Error, ContinueTime does not match header of dispersal file for species %d\n", k);
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
				
				if (col_num < lattice_size)
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


void Ecolattice::loadCompetition(void) {
	/* load the growth and fecundity competition matrices as well as parameters, including species occupancy,
	survival, and dispersal, from the previous simulation. written specifically to work with the way the 
	competition file output is formatted. */

	int col_num;
	int line_num = 0;
	std::string line;
	std::string value;
	std::string name;

	const char *names[] = {"SpeciesOccupancy", "JuvenileSurvival",  "AdultSurvival",  "MaximumCompetition",
						     "DispersalProbability", "DispersalLength", "Fecundity"};

	std::vector<std::vector<double>> parameters(7);
	for (int i = 0; i < 7; i++)
		parameters[i].resize(num_species);

	std::ifstream competition_file;
	competition_file.open(competition_filename);

	if (!competition_file.is_open()) {
			if (id == 0)
				fprintf(stderr, "Error, could not open competition file to load\n");
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
			fprintf(stderr, "Error, competition file ends before anything is loaded\n");
		exit(0);
	}

	if (!getline(competition_file, line)) {
		if (id == 0)
			fprintf(stderr, "Error, competition file ends before anything is loaded\n");
		exit(0);
	}
	line_num++;


	while (getline(competition_file, line)) {
		line_num++;
		col_num = 0;
		if (trimString(line).size() != 0) {
				if (id == 0)
					fprintf(stderr, "Error, formatting issue, line %d in competition file is not a comment\n", line_num);
				exit(0);
		}

		getline(competition_file, line);
		line_num++;

		std::istringstream values(line);

		while (values >> value) {

			col_num++;
			if (col_num > num_species) {
				if (id == 0)
					fprintf(stderr, "Error, too many columns in competition file, line %d\n", line_num);
				exit(0);
			}

			if (col_num == 1) {
				if (value.compare("reset") == 0 && line_num < 17) {
					std::string name(names[ (line_num - 4) / 2 ]);
					int type = 2;
					if (line_num == 4)
						type = 3;
					if (line_num > 12)
						type = 4;
					initializeRandomParameter(parameters[(line_num - 4) / 2], name, type);
					break;
				}
			}

			if (col_num < num_species)
				value = value.substr(0, value.size() - 1);
			
			if (value[0] == '-') {
				fprintf(stderr, "Error, negative value given for competition file line %d column %d\n", line_num, col_num);
				exit(0);
			}

			if (value.find_first_not_of("0123456789.e-") == std::string::npos) {

				try {
					parameters[(line_num - 4) / 2][col_num - 1] = stod(value);
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

		if (line_num == 16)
			break;

	}

	if (line_num < 16) {
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

		std::istringstream values(line);

			while (values >> value) {

				col_num++;
				if (col_num > num_species) {
					if (id == 0)
						fprintf(stderr, "Error, too many columns in competition file, line %d\n", line_num);
					exit(0);
				}			

				if (col_num < num_species)
					value = value.substr(0, value.size() - 1);

				if (value.find_first_not_of("-0123456789.e") == std::string::npos) {
					try {
						competition_fecundity[line_num - 18][col_num - 1] = stod(value);
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
		if (line_num == 17 + num_species)
			break;
	}

	if (line_num < 17 + num_species) {
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

		std::istringstream values(line);

			while (values >> value) {
				col_num++;
				if (col_num > num_species) {
					if (id == 0)
						fprintf(stderr, "Error, too many columns in competition file, line %d\n", line_num);
					exit(0);
				}

				if(col_num < num_species)
					value = value.substr(0, value.size() - 1);

				if (value.find_first_not_of("-0123456789.e") == std::string::npos) {

						try {
							competition_growth[line_num - 19 - num_species][col_num - 1] = stod(value);
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
		if (line_num == 18 + 2 * num_species)
			break;
	}
	if (line_num < 18 + 2 * num_species) {
		if (id == 0)
			fprintf(stderr, "Error, not enough lines from growth competition in competition file, given the number of species\n");
		exit(0);

	}
    for(int i = 0; i < num_species; i++ ) {
		species_occupancy[i] = parameters[0][i];
		juvenile_survival_probability[i] = parameters[1][i];
		adult_survival_probability[i] = parameters[2][i];
		maximum_competition[i] = parameters[3][i];
		dispersal_probability[i] = parameters[4][i];
		dispersal_length[i] = parameters[5][i];
		intrinsic_fecundity[i] = parameters[6][i];
    }
	competition_file.close();
	return;
}

