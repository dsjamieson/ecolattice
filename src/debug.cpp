
	 /**********************************************************
	 * ecolattice
	 *						D.S Jamieson and N.L Kinlock, 2020			
	 *
	 *	driver of the Ecolattice object.  
	 *
	 ***********************************************************/

#include "ecolattice.hpp"

void Ecolattice::compareSpeciesAbundance(void) {
	std::vector<int> counts(num_species, 0);	
	for (int i = 0; i < lattice_size; i++) {
		for (int j = 0; j < lattice_size; j++) {
			if (lattice[i][j]) {
				counts[abs(lattice[i][j]) - 1]++;
			}
		}
	}

	for (int i = 0; i < num_species; i++) {
		fprintf(stdout, "species %d species_abundance = %d\tcounts = %d\n", i, species_abundance[i], counts[i]);
	}
	return;
}


