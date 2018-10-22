
	 /**********************************************************
	 *
	 *			D.S Jamieson and N.L Kinlock, 2018			
	 *
	 *	These methods for the Simulation class read parameters 
	 *	and get simulation variables for ecolattice simultions
	 *
	 ***********************************************************/

#include "simulation.h"

void Simulation::checkInputFormat() {

	const char *names[] = { "Seeds", "BoxSize", "Subdivision", "Species", "Delta", "MaxTimeStep", "InitialOccupancy",
						  "GerminationProbability", "OutfileBase", "OutfileDir", "SpeciesOccupancy", 
						  "JuvenileSurvival", "JuvenileSurvivalSdev", "AdultSurvival", "AdultSurvivalSdev", "MaximumCompetition",
						  "MaximumCompetitionSdev", "DispersalProbability", "DispersalProbabilitySdev", "DispersalLength",
						  "DispersalLengthSdev", "Fecundity", "FecunditySdev", "CompetitionLower", "CompetitionUpper",
						  "CompetitionType", "CompetitionMean", "CompetitionSdev", "CompetitionCorr", "Imbalance", "FecundityTransitivity",
						  "GrowthTransitivity", "RelativeHierarchy", "RestartTime", "CompetitionFile", NULL };

	int i = 0;
	int line_num = 0;
	int pre_trim_size = 0;
	int valid;
	std::string holdstr;
	std::string checkstr;
	std::string line;
	std::ifstream pfile;
	pfile.open(parameter_filename);
	if( pfile.is_open()) {
		while( getline(pfile, line) ){
			line_num++;
			if( line.find("=") != std::string::npos && trimString(line).size() != 0 ) {	
				checkstr = trimStringNoComment(line.substr(0,line.find("=")));
				if( checkstr.find("#") != std::string::npos ) {
					if(id==0)
						fprintf(stderr, "Error, invalid format, # (comment) before = in line %d\n", line_num);
					MPI_Finalize();
					exit(0);
				}
			
				std::istringstream checkstrings( checkstr );
				checkstrings >> holdstr;
				if( checkstrings >> holdstr ) {		
					if(id==0)	
						fprintf(stderr, "Error, invalid format, multiple strings before = in line %d\n", line_num);
					MPI_Finalize();
					exit(0);
				}
				checkstr = trimString(checkstr);
				if( checkstr.size() != 0) {
					valid = 0;
					i = 0;
					while(names[i]!=NULL) {
						if( checkstr.compare(names[i]) == 0 ) {
							valid = 1;
							break;
						}
						i++;
					}
					if(valid==0) {
						if(id==0)
							fprintf(stderr, "Error, unrecognized parameter name %s in line %d\n", checkstr.c_str(), line_num);
						MPI_Finalize();
						exit(0);
					}
				}

			} else if( trimString(line).size() != 0  ) {
				if(id==0)
					fprintf(stderr, "Error, invalid format, line %d does not begin with # (comment) and  no parameter is defined there\n", line_num);
				MPI_Finalize();
				exit(0);
			}
		}
	}
	else{
		if(id==0)
			fprintf(stderr, "Parameter file not found\n");
		MPI_Finalize();
		exit(0);
	}

	return;

}


void Simulation::getSeeds() {

	int i;
	int count = -1;
	int set = 0;

	std::string parameter_name = "Seeds";
	std::string line;
	std::ifstream pfile;
	std::string compstr;
	std::string pstring;
	pfile.open(parameter_filename);
	if( pfile.is_open()) {
		while( getline(pfile, line) ){
			compstr = trimStringNoComment(line.substr(0,line.find("=")));
			if( parameter_name.compare(compstr) == 0 ) {
				set++;
				if(set > 1) {
					if(id==0)
						fprintf(stderr, "Error, multiple lines given for parameter Seeds\n");
					MPI_Finalize();
					exit(0);
				}
				count = 0;
				line = trimString(line.substr( line.find("=")+1, line.length() ));
				std::istringstream pstrings(line);
				while( pstrings >> pstring ) {
					count++;
					pstrings >> std::ws;
				}

				if(count!=0) {
					i=0;
					seeds = new int[count];
					if(!seeds) {
						fprintf(stderr, "Error, unable to allocate memory for seeds\n");
						MPI_Finalize();
						exit(-1);
					}
					pstrings.clear();
					pstrings.seekg(0,pstrings.beg);
					while( pstrings >> pstring ) {
						if(pstring.find_first_not_of( "0123456789" ) == std::string::npos) {
							seeds[i] = stoi(pstring);
							i++;
						}
						else {
							if(id==0)
								fprintf(stderr, "Error, all seeds must be a positive integer\n");
							MPI_Finalize();
							exit(0);
						}
						
						pstrings >> std::ws;
			
					}
				}				
			}
		}
	}
	else{
		if(id==0)
			fprintf(stderr, "Parameter file not found\n");
		MPI_Finalize();
		exit(0);
	}

	if(count==-1){
		if(id==0)
			fprintf(stderr, "Error, found no parameter Seeds\n");
		MPI_Finalize();
		exit(0);
	}


	if(count==0){
		if(id==0)
			fprintf(stderr, "Error, no Seeds given, cannot seed random generator\n");
		MPI_Finalize();
		exit(0);
	}

	seedGenerator(count);

	return;

}


void Simulation::seedGenerator(int num_seeds) {
    std::seed_seq seq(seeds, seeds+num_seeds);
	std::vector<std::uint32_t> seed_vector(std::mt19937::state_size);
    seq.generate(seed_vector.begin(), seed_vector.end());
	std::seed_seq seq2( seed_vector.begin(), seed_vector.end() );
	global_random_generator.seed( seq2 );
	return;
}


void Simulation::getParameter(int *value, std::string parameter_name, int essential ) {

	// Sets value from parameter name in parameter file
	int i;
	int set = 0;
	int check = 0;
	std::string compstr;
	std::string line;
	std::string value_string;
	std::ifstream pfile;
	pfile.open(parameter_filename);

	if( pfile.is_open()) {
		while( getline(pfile, line) ){
			compstr = trimStringNoComment(line.substr(0,line.find("=")));
			if( parameter_name.compare(compstr) == 0 ) {
				check++;
				if(check > 1) {
					if(id==0)
						fprintf(stderr, "Error, multiple lines given for parameter %s\n", parameter_name.c_str());
					MPI_Finalize();
					exit(0);
				}
				line = trimString(line.substr(line.find("=")+1, line.length()));
				std::istringstream pstrings(line);
				if(pstrings >> value_string ) {
					if(value_string.find_first_not_of( "0123456789" ) == std::string::npos) {
						*value = stoi( value_string ) ;
						set = 1;
					}
					else {
						if(id==0)
							fprintf(stderr, "Error, %s must be a positive integer\n", parameter_name.c_str());
						MPI_Finalize();
						exit(0);
					}

					if( pstrings >> value_string) {
						if(id==0)
							fprintf(stderr, "Error, multiple values given for parameter %s\n", parameter_name.c_str());
						MPI_Finalize();
						exit(0);
					}

				}
				else {
					if(essential == 1) {
						if(id==0)
							fprintf(stderr, "Error, no value for %s given\n", parameter_name.c_str());
						MPI_Finalize();
						exit(0);
					}
				}
			}
		}
	}
	else{
		if(id==0)
			fprintf(stderr, "Parameter file not found\n");
		MPI_Finalize();
		exit(0);
	}

	if(set==0 && essential ==1){
		if(id==0)
			fprintf(stderr, "Found no parameter %s\n", parameter_name.c_str());
		MPI_Finalize();
		exit(0);
	}

	return;

}


void Simulation::getParameter(double *value, std::string parameter_name, int essential) {

	// Sets value from parameter name in parameter file
	int i;
	int set = 0;
	int check = 0;
	std::string compstr;
	std::string line;
	std::string value_string;
	std::ifstream pfile;
	pfile.open(parameter_filename);

	if( pfile.is_open()) {
		while( getline(pfile, line) ){
			compstr = trimStringNoComment(line.substr(0,line.find("=")));
			if( parameter_name.compare(compstr) == 0 ) {
				check++;
				if(check > 1) {
					if(id==0)
						fprintf(stderr, "Error, multiple lines given for parameter %s\n", parameter_name.c_str());
					MPI_Finalize();
					exit(0);
				}
				line = trimString(line.substr(line.find("=")+1, line.length()));
				std::istringstream pstrings(line);
				if(pstrings >> value_string) {
					if(value_string.find_first_not_of( "-0123456789.e" ) == std::string::npos) {
						try{
							*value = stod( value_string );
							set = 1;
						}
						catch (...) {
							if(id==0)
								fprintf(stderr, "Error, could not convert value given for parameter %s to double\n", parameter_name.c_str());
							MPI_Finalize();
							exit(0);
						}
					}
					else {
						if(id==0)
							fprintf(stderr, "Error, %s must be a positive integer\n", parameter_name.c_str());
						MPI_Finalize();
						exit(0);
					}

					if( pstrings >> value_string) {
						if(id==0)
							fprintf(stderr, "Error, multiple values given for parameter %s\n", parameter_name.c_str());
						MPI_Finalize();
						exit(0);
					}

				}
				else {
					if(essential == 1) {
						if(id==0)
							fprintf(stderr, "Error, no value for %s given\n", parameter_name.c_str());
						MPI_Finalize();
						exit(0);
					}
				}
			}
		}
	}
	else{
		if(id==0)
			fprintf(stderr, "Parameter file not found\n");
		MPI_Finalize();
		exit(0);
	}

	if(set==0 && essential ==1){
		if(id==0)
			fprintf(stderr, "Found no parameter %s\n", parameter_name.c_str());
		MPI_Finalize();
		exit(0);
	}


	return;

}


void Simulation::getParameter(std::string *value, std::string parameter_name, int essential) {

	// Sets value from parameter name in parameter file

	int i;
	int set = 0;
	int check = 0;
	std::string compstr;
	std::string line;
	std::ifstream pfile;
	pfile.open(parameter_filename);
	if( pfile.is_open()) {
		while( getline(pfile, line) ){
			compstr = trimStringNoComment(line.substr(0,line.find("=")));
			if( parameter_name.compare(compstr) == 0 ) {
				check++;
				if(check > 1) {
					if(id==0)
						fprintf(stderr, "Error, multiple lines given for parameter %s\n", parameter_name.c_str());
					MPI_Finalize();
					exit(0);
				}
				*value = trimString( line.substr(line.find("=")+1, line.length()) ) ;
				if( value->size() == 0 && essential == 1 ) {
					if(id==0)
						fprintf(stderr, "Error, no value for %s given\n", parameter_name.c_str());
					MPI_Finalize();
					exit(0);
				}
				std::istringstream checkstrings( value->c_str() );
				checkstrings >> compstr;
				if( checkstrings >> compstr ) {		
					if(id==0)
						fprintf(stderr, "Error, multiple values given for parameter %s\n", parameter_name.c_str() );
					MPI_Finalize();
					exit(0);
				}
				set = 1;
			}
		}
	}
	else{
		if(id==0)
			fprintf(stderr, "Parameter file not found\n");
		MPI_Finalize();
		exit(0);
	}

	if(set==0 && essential==1){
		if(id==0)
			fprintf(stderr, "Found no parameter %s\n", parameter_name.c_str());
		MPI_Finalize();
		exit(0);
	}

	return;

}


void Simulation::getParameter(double *value_array, int n, std::string parameter_name, int essential) {

	// Sets array values from parameter name in parameter file

	int i;
	int set = 0;
	int check = 0;

	std::string compstr;
	std::string line;
	std::string pstring;
	std::ifstream pfile;
	pfile.open(parameter_filename);
	if( pfile.is_open()) {
		while( getline(pfile, line) ){
			compstr = trimStringNoComment(line.substr(0,line.find("=")));
			if( parameter_name.compare(compstr) == 0 ) {
				check++;
				if(check > 1) {
					if(id==0)
						fprintf(stderr, "Error, multiple lines given for parameter %s\n", parameter_name.c_str());
					MPI_Finalize();
					exit(0);
				}
				i=0;
				line = trimString(line.substr(line.find("=")+1, line.length()));
				std::istringstream pstrings( line );
				while( pstrings >> pstring ) {

					if(i==n+1){
						if(id==0)
							fprintf(stderr, "Error, Number values for %s must be 1 or number of species\n", parameter_name.c_str());
						MPI_Finalize();
						exit(0);
					}

					if(pstring.find_first_not_of( "-0123456789.e" ) == std::string::npos) {

						try{
							value_array[i] = stod(pstring );
						}
						catch (...) {
							if(id==0)
								fprintf(stderr, "Error, could not convert value given for parameter %s to double\n", parameter_name.c_str());
							MPI_Finalize();
							exit(0);
						}

					}
					else {
						if(id==0)
							fprintf(stderr, "Error, %s values must all be floats\n", parameter_name.c_str());
						MPI_Finalize();
						exit(0);					
					}

					if ( ( value_array[i] > 1 || value_array[i] < 0 ) && essential == 2 ) {
						if(id==0)
							fprintf(stderr, "Error, species specific parameter %s treated as probability, all values must be between 0 and 1\n", parameter_name.c_str());
						MPI_Finalize();
						exit(0);
					}
					if ( value_array[i] < 0 && essential == 2 ) {
						if(id==0)
							fprintf(stderr, "Error, species specific parameter %s treated as weight, all values must be non-negative\n", parameter_name.c_str());
						MPI_Finalize();
						exit(0);
					}

					i++;

				} 

				if( i == 0 && essential > 0 ) {
					if(id==0)
						fprintf(stderr, "Error, no value for parameter %s given\n", parameter_name.c_str());
					MPI_Finalize();
					exit(0);
				}			

				if(i==1){
					if(  value_array[0] != 1  && essential == 3 ) {
						if(id==0)
							fprintf(stderr, "Error, species specific parameter %s treated as weights, if only one value is given, it must be 1\n", parameter_name.c_str());
						MPI_Finalize();
						exit(0);
					}
					if ( ( value_array[0] > 1 || value_array[0] < 0 ) && essential == 2 ) {
						if(id==0)
							fprintf(stderr, "Error, species specific parameter %s treated as probability, it must be between 0 and 1\n", parameter_name.c_str());
						MPI_Finalize();
						exit(0);
					}
					for(i=1; i<n; i++)
						value_array[i] = value_array[0];
				}				
				else if( i != n && i!=0  ){
					if(id==0)
						fprintf(stderr, "Error, Number values for %s must be 1 or Species\n", parameter_name.c_str());
					MPI_Finalize();
					exit(0);
				}

				set = 1;

			}
		}
	}
	else{
		if(id==0)
			fprintf(stderr, "Parameter file not found\n");
		MPI_Finalize();
		exit(0);
	}

	if( set==0 && essential > 0 ){
		if(id==0)
			fprintf(stderr, "Found no parameter %s\n", parameter_name.c_str());
		MPI_Finalize();
		exit(0);
	}

	return;

}


void Simulation::initializeNormalRandomArray(double *array, double *mean, double *sdev, int length) {

	// Sets array elements pulled from normal distribtions defined by mean and sdev arrays

	int i;
	for(i=0;i<length;i++){
		std::normal_distribution<float> ndist(mean[i], sdev[i]);
		array[i] = fabs(ndist(generateRandom()));
	}

	return;
}


void Simulation::setRandomParameter(double* parameter_value, int num_species, std::string parameter_name){

	int i;
	double *mean = new double[num_species];
	double *sdev = new double[num_species];
	if(!mean || !sdev) {
		fprintf(stderr, "Error, unable to allocate memory for setting parameter %s\n", parameter_name.c_str());
		MPI_Finalize();
		exit(-1);
	}
	for(i=0; i<num_species; i++) {
		mean[i] = 0.;
		sdev[i] = 0.;
	}

	getParameter(mean, num_species, parameter_name, 1);
	getParameter(sdev, num_species, parameter_name+"Sdev", 0);

	initializeNormalRandomArray(parameter_value, mean, sdev, num_species);

	delete[] mean;
	delete[] sdev;

	return;

}


void Simulation::setRandomProbability(double* parameter_value, int num_species, std::string parameter_name, int type){

	int i;
	double sum = 0.;
	double *mean = new double[num_species];
	double *sdev = new double[num_species];
	if(!mean || !sdev) {
		fprintf(stderr, "Error, unable to allocate memory for setting parameter %s\n", parameter_name.c_str());
		MPI_Finalize();
		exit(-1);
	}
	for(i=0; i<num_species; i++) {
		mean[i] = 0.;
		sdev[i] = 0.;
	}

	getParameter(mean, num_species, parameter_name, type);
	if( type != 3)
		getParameter(sdev, num_species, parameter_name+"Sdev", 0);

	initializeNormalRandomArray(parameter_value, mean, sdev, num_species);

	// type = 3 treats as weights
	if(type == 3) {
		for(i=0; i<num_species; i++)
			sum+=parameter_value[i];

		for(i=0; i<num_species; i++)
			parameter_value[i]/=sum;
	}

	delete[] mean;
	delete[] sdev;

	return;

}


std::string Simulation::trimString(std::string str){

	std::string hold;
	int str_start;
	int str_end;
	int comment_entry;

	str_end = str.size();
	for(comment_entry = 0; comment_entry < str_end; comment_entry++ ){
		if( str[comment_entry]=='#'  )
			break;
	}
	hold = str.substr(0, comment_entry);
	
	str_start = 0;
	str_end = hold.size();

	while(isspace(hold[str_start]))
		str_start++;
	while(isspace(hold[str_end-1]))
		str_end--;

	return hold.substr(str_start, str_end - str_start);

}


std::string Simulation::trimStringNoComment(std::string str){

	int str_start;
	int str_end;

	str_start = 0;
	str_end = str.size();

	while(isspace(str[str_start]))
		str_start++;
	while(isspace(str[str_end-1]))
		str_end--;

	return str.substr(str_start, str_end - str_start);

} 

