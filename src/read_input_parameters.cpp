
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2020			
	 *
	 *		methods for the Simulation class. read in parameters 
	 *		used in simulations from file.
	 *
	 ***********************************************************/

#include "ecolattice.hpp"

void Ecolattice::checkInputFormat(void) {
	/* check 'parameters.dat' for all parameters and correct format. can remove commented lines from parameters file
		where comments are specified with '#'. */

	const char *names[] = {"LatticeSize", "Subdivision", "Species", "Delta", "MaxTimeStep", "InitialOccupancy",
						  "GerminationProbability", "OutfileBase", "OutfileDir", "SpeciesOccupancy", "SpeciesOccupancySdev",
						  "JuvenileSurvival", "JuvenileSurvivalSdev", "AdultSurvival", "AdultSurvivalSdev", "MaximumCompetition",
						  "MaximumCompetitionSdev", "DispersalProbability", "DispersalProbabilitySdev", "DispersalLength",
						  "DispersalLengthSdev", "Fecundity", "FecunditySdev", "CompetitionLower", "CompetitionUpper",
						  "CompetitionDiagLower", "CompetitionDiagUpper", "CompetitionDiagMean", "CompetitionDiagSdev",
						  "CompetitionType", "CompetitionMean", "CompetitionSdev", "CompetitionCorr", "Imbalance", "FecundityTransitivity",
						  "GrowthTransitivity", "RelativeHierarchy", "ContinueTime", "CompetitionFile", "MinPersistence", "NumThreads", NULL};

	int i = 0;
	int line_num = 0;
	int valid;
	std::string holdstr;
	std::string checkstr;
	std::string line;
	std::ifstream pfile;
	pfile.open(parameter_filename);

	if (!pfile.is_open()) {
		if (id == 0)
			fprintf(stderr, "Parameter file not found\n");
		exit(0);
	}

	while (getline(pfile, line)) {
		line_num++;
		// comments allowed only after parameters are specified
		if (line.find("=") != std::string::npos && trimString(line).size() != 0) {	
			checkstr = trimStringNoComment(line.substr(0, line.find("=")));
			if (checkstr.find("#") != std::string::npos) {
				if (id == 0)
					fprintf(stderr, "Error, invalid format, # (comment) before = in line %d\n", line_num);
				exit(0);
			}
		
			std::istringstream checkstrings(checkstr);
			checkstrings >> holdstr;
			// only one name for each parameter, no spaces allowed in parameter names
			if (checkstrings >> holdstr) {		
				if (id == 0)	
					fprintf(stderr, "Error, invalid format, multiple strings before = in line %d\n", line_num);
				exit(0);
			}
			checkstr = trimString(checkstr);
			// check parameter names against allowed names
			if (checkstr.size() != 0) {
				valid = 0;
				i = 0;
				while (names[i] != NULL) {
					if (checkstr.compare(names[i]) == 0) {
						valid = 1;
						break;
					}
					i++;
				}
				if (valid == 0) {
					if (id == 0)
						fprintf(stderr, "Error, unrecognized parameter name %s in line %d\n", checkstr.c_str(), line_num);
					exit(0);
				}
			}
		} else if (trimString(line).size() != 0) {
			if (id == 0)
				fprintf(stderr, "Error, invalid format in  line %d. Lines must begin with # (comment) or a parameter name\n", line_num);
			exit(0);
		}
	}

	pfile.close();

	return;
}

int Ecolattice::getParameter(int & value, std::string parameter_name, int essential ) {
	/* sets value of parameter given the parameter name in parameter file, 'parameter.dat' 
	file is already open from previous method, checkInputFormat. method for integer parameters. */

	int set = 0;
	int check = 0;
	std::string compstr;
	std::string line;
	std::string value_string;
	std::ifstream pfile;
	pfile.open(parameter_filename);

	if (!pfile.is_open()) {
		if (id == 0)
			fprintf(stderr, "Parameter file not found\n");
		exit(0);
	}		
	
	while (getline(pfile, line) ){
		compstr = trimStringNoComment(line.substr(0, line.find("=")));
		if (parameter_name.compare(compstr) == 0) {
			check++;
			if (check > 1) {
				if (id == 0)
					fprintf(stderr, "Error, multiple lines given for parameter %s\n", parameter_name.c_str());
				exit(0);
			}
			line = trimString(line.substr(line.find("=") + 1, line.length()));
			std::istringstream pstrings(line);
			if (pstrings >> value_string) {
				if (value_string.find_first_not_of("0123456789") == std::string::npos) {
					value = stoi(value_string);
					set = 1;
				}
				else {
					if (id == 0)
						fprintf(stderr, "Error, %s must be a positive integer\n", parameter_name.c_str());
					exit(0);
				}
				if (pstrings >> value_string) {
					if (id == 0)
						fprintf(stderr, "Error, multiple values given for parameter %s\n", parameter_name.c_str());
					exit(0);
				}
			}
			else {
				if (essential == 1) {
					if (id == 0)
						fprintf(stderr, "Error, no value for %s given\n", parameter_name.c_str());
					exit(0);
				}
			}
		}
	}

	if (set == 0 && essential == 1){
		if (id == 0)
			fprintf(stderr, "Found no parameter %s\n", parameter_name.c_str());
		exit(0);
	}
	pfile.close();
	return set;
}

void Ecolattice::getParameter(double & value, std::string parameter_name, int essential) {
	/* sets value of parameter given the parameter name in parameter file, 'parameter.dat' 
	file is already open from previous method, checkInputFormat. method for double parameters. */

	int set = 0;
	int check = 0;
	std::string compstr;
	std::string line;
	std::string value_string;
	std::ifstream pfile;
	pfile.open(parameter_filename);

	if (!pfile.is_open()) {
		if (id == 0)
			fprintf(stderr, "Parameter file not found\n");
		exit(0);
	}
	
	while (getline(pfile, line)) {
		compstr = trimStringNoComment(line.substr(0, line.find("=")));
		if (parameter_name.compare(compstr) == 0) {
			check++;
			if (check > 1) {
				if (id == 0)
					fprintf(stderr, "Error, multiple lines given for parameter %s\n", parameter_name.c_str());
				exit(0);
			}
			line = trimString(line.substr(line.find("=") + 1, line.length()));
			std::istringstream pstrings(line);
			if (pstrings >> value_string) {
				if (value_string.find_first_not_of("-0123456789.e") == std::string::npos) {
					try {
						value = stod(value_string);
						set = 1;
					}
					catch (...) {
						if (id == 0)
							fprintf(stderr, "Error, could not convert value given for parameter %s to double\n", parameter_name.c_str());
						exit(0);
					}
				}
				else {
					if (id == 0)
						fprintf(stderr, "Error, %s must be a double\n", parameter_name.c_str());
					exit(0);
				}
				if (pstrings >> value_string) {
					if (id == 0)
						fprintf(stderr, "Error, multiple values given for parameter %s\n", parameter_name.c_str());
					exit(0);
				}

			}
			else {
				if (essential == 1) {
					if (id == 0)
						fprintf(stderr, "Error, no value for %s given\n", parameter_name.c_str());
					exit(0);
				}
			}
		}
	}

	if (set == 0 && essential == 1){
		if (id == 0)
			fprintf(stderr, "Found no parameter %s\n", parameter_name.c_str());
		exit(0);
	}

	pfile.close();
	return;

}

void Ecolattice::getParameter(std::string & value, std::string parameter_name, int essential) {
	/* sets value of parameter given the parameter name in parameter file, 'parameter.dat' 
	file is already open from previous method, checkInputFormat. method for string parameters. */

	int set = 0;
	int check = 0;
	std::string compstr;
	std::string line;
	std::ifstream pfile;
	pfile.open(parameter_filename);

	if (!pfile.is_open()) {
		if (id == 0)
			fprintf(stderr, "Parameter file not found\n");
		exit(0);
	}

	while (getline(pfile, line)){
		compstr = trimStringNoComment(line.substr(0, line.find("=")));
		if (parameter_name.compare(compstr) == 0) {
			check++;
			if (check > 1) {
				if (id == 0)
					fprintf(stderr, "Error, multiple lines given for parameter %s\n", parameter_name.c_str());
				exit(0);
			}
			value = trimString(line.substr(line.find("=") + 1, line.length())) ;
			if (value.size() == 0 && essential == 1) {
				if (id==0)
					fprintf(stderr, "Error, no value for %s given\n", parameter_name.c_str());
				exit(0);
			}
			std::istringstream checkstrings(value.c_str());
			checkstrings >> compstr;
			if (checkstrings >> compstr) {		
				if (id == 0)
					fprintf(stderr, "Error, multiple values given for parameter %s\n", parameter_name.c_str() );
				exit(0);
			}
			set = 1;
		}
	}

	if (set == 0 && essential == 1) {
		if (id == 0)
			fprintf(stderr, "Found no parameter %s\n", parameter_name.c_str());
		exit(0);
	}
	pfile.close();
	return;

}


void Ecolattice::getParameter(std::vector<int> & value_array, std::string parameter_name, int essential) {
	/* sets value of parameter given the name in parameter file, which is an argument for Simulation.
	file is already open from previous method, checkInputFormat. method for array parameters.
	this method has some unique requirements, some arrays require positive integers [0, ) (essential
	type 2) and will return an error if not satisfied. */

	unsigned long i;
	int set = 0;
	int check = 0;

	std::string compstr;
	std::string line;
	std::string pstring;
	std::ifstream pfile;
	pfile.open(parameter_filename);

	if (!pfile.is_open()) {
		if (id == 0)
			fprintf(stderr, "Parameter file not found\n");
		exit(0);
	}

	while (getline(pfile, line)){
		compstr = trimStringNoComment(line.substr(0, line.find("=")));
		if (parameter_name.compare(compstr) == 0) {
			check++;
			if (check > 1) {
				if (id == 0)
					fprintf(stderr, "Error, multiple lines given for parameter %s\n", parameter_name.c_str());
				exit(0);
			}
			i = 0;
			line = trimString(line.substr(line.find("=") + 1, line.length()));
			std::istringstream pstrings(line);
			while (pstrings >> pstring) {
				if (i == value_array.size() + 1){
					if (id == 0)
						fprintf(stderr, "Error, Number of values for parameter %s must be 1 or number of species\n", parameter_name.c_str());
					exit(0);
				}

				if (pstring.find_first_not_of("-0123456789") == std::string::npos) {

					try {
						value_array[i] = stoi(pstring);
					}
					catch (...) {
						if (id == 0)
							fprintf(stderr, "Error, could not convert value given for parameter %s to int\n", parameter_name.c_str());
						exit(0);
					}
				}
				else {
					if (id == 0)
						fprintf(stderr, "Error, %s values must all be integers\n", parameter_name.c_str());
					exit(0);					
				}

				if (value_array[i] < 0 && essential == 2) {
					if (id == 0)
						fprintf(stderr, "Error, species specific parameter %s must be non-negative integer\n", parameter_name.c_str());
					exit(0);
				}
			

				i++;

			}
 
			if (i == 0 && essential > 0) {
				if (id == 0)
					fprintf(stderr, "Error, no value for parameter %s given\n", parameter_name.c_str());
				exit(0);
			}			
			if (i == 1){

				for (i = 1; i < value_array.size(); i++)
					value_array[i] = value_array[0];

			}				
			else if (i != value_array.size() && i!= 0){
				if (id == 0)
					fprintf(stderr, "Error, Number values for %s must be 1 or Species\n", parameter_name.c_str());
				exit(0);
			}
			set = 1;
		}
	}	

	if (set == 0 && essential > 0) {
		if (id == 0)
			fprintf(stderr, "Found no parameter %s\n", parameter_name.c_str());
		exit(0);
	}
	pfile.close();
	return;
}

void Ecolattice::getParameter(std::vector<double> & value_array, std::string parameter_name, int essential) {
	/* sets value of parameter given the name in parameter file, which is an argument for Simulation.
	file is already open from previous method, checkInputFormat. method for array parameters.
	this method has some unique requirements, some arrays require probabilities [0, 1] (essential
	type 2), weights [0, ) (type 3), or non-negative values [0, ) (type 4) , and will return errors if not satisfied. */

	unsigned long i;
	int set = 0;
	int check = 0;
	double value_sum = 0.;

	std::string compstr;
	std::string line;
	std::string pstring;
	std::ifstream pfile;
	pfile.open(parameter_filename);

	if (!pfile.is_open()) {
		if (id == 0)
			fprintf(stderr, "Parameter file not found\n");
		exit(0);
	}

	while (getline(pfile, line)){
		compstr = trimStringNoComment(line.substr(0, line.find("=")));
		if (parameter_name.compare(compstr) == 0) {
			check++;
			if (check > 1) {
				if (id == 0)
					fprintf(stderr, "Error, multiple lines given for parameter %s\n", parameter_name.c_str());
				exit(0);
			}
			i = 0;
			line = trimString(line.substr(line.find("=") + 1, line.length()));
			std::istringstream pstrings(line);
			while (pstrings >> pstring) {
				if (i == value_array.size() + 1){
					if (id == 0)
						fprintf(stderr, "Error, Number values for %s must be 1 or number of species\n", parameter_name.c_str());
					exit(0);
				}

				if (pstring.find_first_not_of("-0123456789.e") == std::string::npos) {

					try {
						value_array[i] = stod(pstring);
						value_sum += value_array[i];
					}
					catch (...) {
						if (id == 0)
							fprintf(stderr, "Error, could not convert value given for parameter %s to double\n", parameter_name.c_str());
						exit(0);
					}
				}
				else {
					if (id == 0)
						fprintf(stderr, "Error, %s values must all be floats\n", parameter_name.c_str());
					exit(0);					
				}

				if ((value_array[i] > 1 || value_array[i] < 0 ) && essential == 2) {
					if (id == 0)
						fprintf(stderr, "Error, species specific parameter %s treated as probability, all values must be between 0 and 1\n", parameter_name.c_str());
					exit(0);
				}
				if (value_array[i] < 0 && essential == 3) {
					if (id == 0)
						fprintf(stderr, "Error, species specific parameter %s treated as weight, all values must be non-negative\n", parameter_name.c_str());
					exit(0);
				}
				if (value_array[i] < 0 && essential == 4) {
					if (id == 0)
						fprintf(stderr, "Error, species specific parameter %s must be non-negative\n", parameter_name.c_str());
					exit(0);
				}

				i++;

			}
 
			if (i == 0 && essential > 0) {
				if (id == 0)
					fprintf(stderr, "Error, no value for parameter %s given\n", parameter_name.c_str());
				exit(0);
			}			
			if (i == 1){
				if (value_array[0] != 1  && essential == 3) {
					if (id == 0)
						fprintf(stderr, "Error, species specific parameter %s treated as weights, if only one value is given, it must be 1\n", parameter_name.c_str());
					exit(0);
				}

				for (i = 1; i < value_array.size(); i++)
					value_array[i] = value_array[0];

				value_sum = value_array[0]*value_array.size();

			}				
			else if (i != value_array.size() && i!= 0){
				if (id == 0)
					fprintf(stderr, "Error, Number values for %s must be 1 or Species\n", parameter_name.c_str());
				exit(0);
			}
			set = 1;
		}
	}

	// values treated as weights if essential type = 3
	if(essential == 3) {
		for(i = 0; i < value_array.size(); i++)
			value_array[i] /= value_sum;
	}
	

	if (set == 0 && essential > 0) {
		if (id == 0)
			fprintf(stderr, "Found no parameter %s\n", parameter_name.c_str());
		exit(0);
	}
	pfile.close();
	return;
}

std::string Ecolattice::trimString(std::string str) {
	std::string hold;
	int str_start, str_end, comment_entry;
	str_end = str.size();
	for (comment_entry = 0; comment_entry < str_end; comment_entry++) {
		if (str[comment_entry]=='#')
			break;
	}
	hold = str.substr(0, comment_entry);
	str_start = 0;
	str_end = hold.size();
	while (isspace(hold[str_start]))
		str_start++;
	while (isspace(hold[str_end - 1]))
		str_end--;
	return hold.substr(str_start, str_end - str_start);
}

std::string Ecolattice::trimStringNoComment(std::string str) {
	int str_start;
	int str_end;
	str_start = 0;
	str_end = str.size();
	while (isspace(str[str_start]))
		str_start++;
	while (isspace(str[str_end - 1]))
		str_end--;
	return str.substr(str_start, str_end - str_start);
} 

