#include "fluid2d.h"
//#include<Eigen/Core>
#include <string>
#include <cmath>
#include<iostream>


// definitions of methods for the class Fluid2d are located
// methods for boundary condition





// reflecting boundary condition of selected type
template<int Ni, int Nj, int Nb>
void Fluid2d<Ni, Nj, Nb>::reflect_bc(const std::string& bc_type) {
	if (bc_type == "periodical_in_i") {
		bc_periodical_in_i();
	}
	else {
		std::cout << "ERROR Fluid2d::reflect_bc > BC-type cannnot be specified.";
		exit(EXIT_FAILURE);
	}
};





//------------------Procedures for Each BC-Type-----------------------

// periodical in i-direction, Dirichlet on j-boundaries
template<int Ni, int Nj, int Nb>
void Fluid2d<Ni, Nj, Nb>::bc_periodical_in_i() {
	// preiodical in i-direction
	// for i=0
	for (int i = 0; i < Nb; ++i) {
		for (int j = 0; j < Nj; ++j) {
			rho[i][j] = rho[Ni-2*Nb+i][j];
			u[i][j] = u[Ni-2*Nb+i][j];
			v[i][j] = v[Ni-2*Nb+i][j];
			e[i][j] = e[Ni-2*Nb+i][j];
		}
	}
	// for i=Ni
	for (int i = 0; i < Nb; ++i) {
		for (int j = 0; j < Nj; ++j) {
			rho[Ni-Nb+i][j] = rho[Nb+i][j];
			u[Ni-Nb+i][j] = u[Nb+i][j];
			v[Ni-Nb+i][j] = v[Nb+i][j];
			e[Ni-Nb+i][j] = e[Nb+i][j];
		}
	}
	// No updates for j-boundaries since Dirichlet
};