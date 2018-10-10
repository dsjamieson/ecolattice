// Percentage juvenile in initial conditions option
// optaional condiiton function for parameter file
// initialization defaults
// work out overfilling box problem
#include <fstream>
#include <sstream>
#include <random>
#include "simulation.h"

void Simulation::getSeeds() {

	int i;
	int count = 0;

	std::string parameter_name = "Seeds";
	std::string line;
	std::ifstream pfile;
	pfile.open(parameter_filename);
	if( pfile.is_open()) {
		while( getline(pfile, line) ){
			if( parameter_name.compare(line.substr(0,line.find(" = "))) == 0 ) {
				std::string pstring;
				std::istringstream pstrings( line.substr(line.find("=")+2, line.length())  );
				while( getline( pstrings, pstring, ' ')  )
					count++;
			
				if(count!=0) {
					i=0;
					seeds = new int[count];
					pstrings.clear();
					pstrings.seekg(0,pstrings.beg);
					while( getline( pstrings, pstring, ' ')  ) {
						seeds[i] = stoi(pstring);
						i++;
					}

				}
				
				break;
			}
		}
	}
	else{
		fprintf(stderr, "Parameter file not found\n");
		exit(0);
	}

	if(count==0){
		fprintf(stderr, "Error, no Seeds given, cannot seed random generator\n");
			exit(0);
	}

	seedGenerator(count);

	return;

}

void Simulation::getParameter(int *value, std::string parameter_name, int essential ) {

	// Sets value from parameter name in parameter file

	int i, set = 0;
	std::string compstr;
	std::string line;
	std::string value_string;
	std::ifstream pfile;
	pfile.open(parameter_filename);
	if( pfile.is_open()) {
		while( getline(pfile, line) ){
			compstr = trimstr(line.substr(0,line.find("=")));
			if( parameter_name.compare(compstr) == 0 ) {
				value_string = trimstr(line.substr(line.find("=")+1, line.length()));
				if(value_string.find_first_not_of( "0123456789" ) == std::string::npos) {
					*value = stoi( value_string ) ;
					set = 1;
					break;
				}
				else {
					fprintf(stderr, "Error, %s must be and integer\n", parameter_name.c_str());
					exit(0);
				}
			}
		}
	}
	else{
		fprintf(stderr, "Parameter file not found\n");
		exit(0);
	}

	if(set==0){
		if(essential) {
			fprintf(stderr, "Found no parameter %s\n", parameter_name.c_str());
			exit(0);
		}
	}

	return;

}


void Simulation::getParameter(double *value, std::string parameter_name, int essential) {

	// Sets value from parameter name in parameter file

	int i, set = 0;
	std::string compstr;
	
	std::string line;
	std::ifstream pfile;
	pfile.open(parameter_filename);
	if( pfile.is_open()) {
		while( getline(pfile, line) ){
			compstr = trimstr(line.substr(0,line.find("=")));
			if( parameter_name.compare(compstr) == 0 ) {
				*value = stod( line.substr(line.find("=")+1, line.length()) ) ;
				set = 1;
				break;
			}
		}
	}
	else{
		if(essential) {
			fprintf(stderr, "Parameter file not found\n");
			exit(0);
		}
	}

	if(set==0){
		fprintf(stderr, "Found no parameter %s\n", parameter_name.c_str());
		exit(0);
	}

	return;

}

void Simulation::getParameter(std::string *value, std::string parameter_name, int essential) {

	// Sets value from parameter name in parameter file

	int i, set = 0;
	std::string compstr;
	
	std::string line;
	std::ifstream pfile;
	pfile.open(parameter_filename);
	if( pfile.is_open()) {
		while( getline(pfile, line) ){
			compstr = trimstr(line.substr(0,line.find("=")));
			if( parameter_name.compare(compstr) == 0 ) {
				*value = trimstr( line.substr(line.find("=")+1, line.length()) ) ;
				set = 1;
				break;
			}
		}
	}
	else{
		if(essential) {
			fprintf(stderr, "Parameter file not found\n");
			exit(0);
		}
	}

	if(set==0){
		fprintf(stderr, "Found no parameter %s\n", parameter_name.c_str());
		exit(0);
	}

	return;

}

void Simulation::getParameter(double *value_array, int n, std::string parameter_name, int essential) {

	// Sets array values from parameter name in parameter file

	int i;
	int set = 0;

	std::string line;
	std::ifstream pfile;
	pfile.open(parameter_filename);
	if( pfile.is_open()) {
		while( getline(pfile, line) ){
			if( parameter_name.compare(line.substr(0,line.find(" = "))) == 0 ) {
				i=0;
				std::string pstring;
				std::istringstream pstrings( line.substr(line.find("=")+2, line.length())  );
				while( getline( pstrings, pstring, ' ')  ){
					if(i==n+1){
						fprintf(stderr, "Error, Number values for %s must be 1 or number of species\n", parameter_name.c_str());
						exit(0);
					}
					value_array[i] = stod( pstring ) ;
					i++;
				}
				if(i==1){
					for(i=1; i<n; i++)
						value_array[i] = value_array[0];
				}				
				if(i<n){
					fprintf(stderr, "Error, Number values for %s must be 1 or number of species\n", parameter_name.c_str());
					exit(0);
				}
				set = 1;
				break;
			}
		}
	}
	else{
		fprintf(stderr, "Parameter file not found\n");
		exit(0);
	}

	if(set==0){
		if(essential) {
			fprintf(stderr, "Found no parameter %s\n", parameter_name.c_str());
			exit(0);
		}
	}

	return;

}

void Simulation::initializeNormalRandomArray(double *array, double *mean, double *sdev, int length) {

	// Sets array elements pulled from normal distribtions defined by mean and sdev arrays

	int i;
	for(i=0;i<length;i++){
		std::normal_distribution<float> ndist(mean[i], sdev[i]);
		array[i] = fabs(ndist(global_random_generator));
	}

	return;
}

void Simulation::setRandomParameter(double* parameter_value, int num_species, std::string parameter_name){

	double * mean = new double[num_species];
	double * sdev = new double[num_species];

	getParameter(mean, num_species, parameter_name, 1);
	getParameter(sdev, num_species, parameter_name+"Sdev", 1);

	initializeNormalRandomArray(parameter_value, mean, sdev, num_species);

	return;

}



std::string Simulation::trimstr(std::string str){

	int start = 0;
	int end = str.size();
	while(isspace(str[start]))
		start++;
	while(isspace(str[end-1]))
		end--;
	
	return str.substr(start, end);

} 
