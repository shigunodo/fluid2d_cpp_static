#include "fluid2d.h"
#include "fnd.h" // declaration of class Fnd is here
#include <cmath>

// definitions of methods for the class Fluid2d are located
// methods for time marching, calculation of RHS, reconstruction scheme




//=============================================================
//------------------time marching schemes----------------------
//=============================================================


//-------------3rd order SSP Runge-Kutta method----------------
// Shu & Osher (1988)
template<int Ni, int Nj, int Nb>
void Fluid2d<Ni, Nj, Nb>::march_ssprk3(const double& dt, 
	const std::string& bc_type,
	const std::string& reconstruction,
	const std::string& flux_scheme) {
	double SA;
	constexpr int Nf = 4; // components of flux/conserved var


	// construction of conserved var
	for (int i = 0; i < Ni-2*Nb; ++i) {
		for (int j = 0; j < Nj-2*Nb; ++j) {
			SA = 0.25 * (S[Nb + i - 1][Nb + j - 1] + S[Nb + i][Nb + j - 1] \
				+ S[Nb + i - 1][Nb + j] + S[Nb + i][Nb + j]);
			calc_conserved(rho[Nb + i][Nb + j], u[Nb + i][Nb + j], \
				v[Nb + i][Nb + j], e[Nb + i][Nb + j], SA, arr_Q0[i][j]);
		}
	}


	// 1st stage
	// calc rhs
	rhs(rho, u, v, e, arr_Q1, reconstruction, flux_scheme);
	// marching
	for (int i = 0;i < Ni-2*Nb; ++i) {
		for (int j = 0; j < Nj-2*Nb; ++j) {
			for (int k = 0; k < Nf; ++k) {
				arr_Q1[i][j][k] = arr_Q0[i][j][k] + dt * arr_Q1[i][j][k];
			}
			// convert into basic var
			SA = 0.25 * (S[Nb + i - 1][Nb + j - 1] + S[Nb + i][Nb + j - 1] \
				+ S[Nb + i - 1][Nb + j] + S[Nb + i][Nb + j]);
			calc_basic(arr_Q1[i][j], SA, rho[Nb + i][Nb + j], 
				u[Nb + i][Nb + j], v[Nb + i][Nb + j], e[Nb + i][Nb + j]);
		}
	}
	// reflecting boundary condition
	reflect_bc(bc_type);


	// 2nd stage
	// calc rhs
	rhs(rho, u, v, e, arr_Q2, reconstruction, flux_scheme);
	// marching
	for (int i = 0;i < Ni-2*Nb; ++i) {
		for (int j = 0; j < Nj-2*Nb; ++j) {
			for (int k = 0; k < Nf; ++k) {
				arr_Q2[i][j][k] = 0.75 * arr_Q0[i][j][k] \
				 + 0.25 * (arr_Q1[i][j][k] + dt * arr_Q2[i][j][k]);
			}
			// convert into basic var
			SA = 0.25 * (S[Nb + i - 1][Nb + j - 1] + S[Nb + i][Nb + j - 1] \
				+ S[Nb + i - 1][Nb + j] + S[Nb + i][Nb + j]);
			calc_basic(arr_Q2[i][j], SA, rho[Nb + i][Nb + j], 
				u[Nb + i][Nb + j], v[Nb + i][Nb + j], e[Nb + i][Nb + j]);
		}
	}
	// reflecting boundary condition
	reflect_bc(bc_type);


	// 3rd stage
	// calc rhs
	rhs(rho, u, v, e, arr_Q1, reconstruction, flux_scheme);
	// marching
	for (int i = 0;i < Ni-2*Nb; ++i) {
		for (int j = 0; j < Nj-2*Nb; ++j) {
			for (int k = 0; k < Nf; ++k) {
				arr_Q0[i][j][k] = arr_Q0[i][j][k] / 3.0 \
				 + 2.0 / 3.0 * (arr_Q2[i][j][k] + dt * arr_Q1[i][j][k]);
			}
			// convert into basic var
			SA = 0.25 * (S[Nb + i - 1][Nb + j - 1] + S[Nb + i][Nb + j - 1] \
				+ S[Nb + i - 1][Nb + j] + S[Nb + i][Nb + j]);
			calc_basic(arr_Q0[i][j], SA, rho[Nb + i][Nb + j], 
				u[Nb + i][Nb + j], v[Nb + i][Nb + j], e[Nb + i][Nb + j]);
		}
	}
	// reflecting boundary condition
	reflect_bc(bc_type);



	// time
	t += dt;
	tstep += 1;
};





//=============================================================
//----------------------CFL Condition--------------------------
//=============================================================
// dt = CFL_coeff * min( dx / (cs + sqrt(u*u+v*v)) )
// where cs is local sound speed
template<int Ni, int Nj, int Nb>
double Fluid2d<Ni, Nj, Nb>::calc_cfl(const double CFL_coeff) const {
	double nu = 1.0e+10;
	double tmp;
	for (int i = 0; i < Ni - 2 * Nb; ++i) {
		for (int j = 0; j < Nj - 2 * Nb; ++j) {
			tmp = dx[i][j] \
				/ (calc_cs(rho[i + Nb][j + Nb], u[i + Nb][j + Nb],
					v[i + Nb][j + Nb], e[i + Nb][j + Nb]) \
					+ std::sqrt(u[i + Nb][j + Nb] * u[i + Nb][j + Nb] \
						+ v[i + Nb][j + Nb] * v[i + Nb][j + Nb]));
			nu = std::min(nu, tmp);
		}
	}
	return CFL_coeff * nu;
};




//=============================================================
//----------------Evaluating RHS of Equation-------------------
//=============================================================
template<int Ni, int Nj, int Nb>
void Fluid2d<Ni, Nj, Nb>::rhs(const double(&rho)[Ni][Nj], 
	const double(&u)[Ni][Nj],
	const double(&v)[Ni][Nj], const double(&e)[Ni][Nj],
	double(&arr_Q)[Ni-2*Nb][Nj-2*Nb][Nf], 
	const std::string& reconstruction,
	const std::string& flux_scheme) {

	// temporal variables
	double SA, ixSA, iySA, jxSA, jySA;
	double rhoL, uL, vL, eL, rhoR, uR, vR, eR;


	//---------------------i-direction-------------------------
	// evaluating numerical flux at cell-boundary (i+0.5,j)

	for (int i = 0; i < Ni-2*Nb+1; ++i) {
		for (int j = 0; j < Nj-2*Nb; ++j) {

			// inverse of Jacobian
			// evaluated at (i + 0.5, j)
			// by averaging the values at (i+0.5,j+0.5) and (i+0.5,j-0.5)
			SA = 0.5 * (S[Nb + i - 1][Nb + j - 1] + S[Nb + i - 1][Nb + j]);

			// metrices devided by Jacobian
			// evaluated at (i+0.5,j)
			// by averaging the values at (i,j) and (i+1,j)
			ixSA = 0.5 * (ixS[i][j + 1] + ixS[i + 1][j + 1]);
			iySA = 0.5 * (iyS[i][j + 1] + iyS[i + 1][j + 1]);
			jxSA = 0.5 * (jxS[i][j + 1] + jxS[i + 1][j + 1]);
			jySA = 0.5 * (jyS[i][j + 1] + jyS[i + 1][j + 1]);

			// reconstruction
			if (reconstruction == "MUSCL_minmod_basic") {
				reconstruction_basic_muscl3( \
					rho[Nb+i-2][Nb+j], rho[Nb+i-1][Nb+j], 
					rho[Nb+i][Nb+j], rho[Nb+i+1][Nb+j],
					u[Nb+i-2][Nb+j], u[Nb+i-1][Nb+j], 
					u[Nb+i][Nb+j], u[Nb+i+1][Nb+j],
					v[Nb+i-2][Nb+j], v[Nb+i-1][Nb+j], 
					v[Nb+i][Nb+j], v[Nb+i+1][Nb+j],
					e[Nb+i-2][Nb+j], e[Nb+i-1][Nb+j], 
					e[Nb+i][Nb+j], e[Nb+i+1][Nb+j],
					rhoL, uL, vL, eL, rhoR, uR, vR, eR);
			}
			else if (reconstruction == "MP5_basic") {
				reconstruction_basic_mp5(\
					rho[Nb+i-3][Nb+j], rho[Nb+i-2][Nb+j],
					rho[Nb+i-1][Nb+j], rho[Nb+i][Nb+j],
					rho[Nb+i+1][Nb+j], rho[Nb+i+2][Nb+j],
					u[Nb+i-3][Nb+j], u[Nb+i-2][Nb+j],
					u[Nb+i-1][Nb+j], u[Nb+i][Nb+j],
					u[Nb+i+1][Nb+j], u[Nb+i+2][Nb+j],
					v[Nb+i-3][Nb+j], v[Nb+i-2][Nb+j],
					v[Nb+i-1][Nb+j], v[Nb+i][Nb+j],
					v[Nb+i+1][Nb+j], v[Nb+i+2][Nb+j],
					e[Nb+i-3][Nb+j], e[Nb+i-2][Nb+j],
					e[Nb+i-1][Nb+j], e[Nb+i][Nb+j],
					e[Nb+i+1][Nb+j], e[Nb+i+2][Nb+j],
					rhoL, uL, vL, eL, rhoR, uR, vR, eR);
			}
			else {
				std::cout << "ERROR Fluid2d::rhs > Reconstruction scheme cannnot be specified.";
				exit(EXIT_FAILURE);
			}

			// evaluating flux using flux scheme
			if (flux_scheme == "Roe_FDS") {
				roe_fds(rhoL, uL, vL, eL, rhoR, uR, vR, eR,
					ixSA, iySA, SA, arr_Fi[i][j]);
			}
			else {
				std::cout << "ERROR Fluid2d::rhs > Flux scheme cannnot be specified.";
				exit(EXIT_FAILURE);
			}
		} // for j
	} // for i


	//---------------------j-direction-------------------------
	// evaluating numerical flux at (i,j+0.5)

	for (int i = 0; i < Ni-2*Nb; ++i) {
		for (int j = 0; j < Nj-2*Nb+1; ++j) {

			// inverse of Jacobian
			// evaluated at(i, j+0.5)
			// by averaging the values at(i+0.5,j+0.5) and (i-0.5,j+0.5)
			SA = 0.5 * (S[Nb + i - 1][Nb + j - 1] + S[Nb + i][Nb + j - 1]);

			// metrices devided by Jacobian
			// evaluated at (i,j+0.5)
			// by averaging the values at (i,j) and (i,j+1)
			ixSA = 0.5 * (ixS[i + 1][j] + ixS[i + 1][j + 1]);
			iySA = 0.5 * (iyS[i + 1][j] + iyS[i + 1][j + 1]);
			jxSA = 0.5 * (jxS[i + 1][j] + jxS[i + 1][j + 1]);
			jySA = 0.5 * (jyS[i + 1][j] + jyS[i + 1][j + 1]);

			// reconstruction
			if (reconstruction == "MUSCL_minmod_basic") {
				reconstruction_basic_muscl3( \
					rho[Nb+i][Nb+j-2], rho[Nb+i][Nb+j-1], 
					rho[Nb+i][Nb+j], rho[Nb+i][Nb+j+1],
					u[Nb+i][Nb+j-2], u[Nb+i][Nb+j-1], 
					u[Nb+i][Nb+j], u[Nb+i][Nb+j+1],
					v[Nb+i][Nb+j-2], v[Nb+i][Nb+j-1], 
					v[Nb+i][Nb+j], v[Nb+i][Nb+j+1],
					e[Nb+i][Nb+j-2], e[Nb+i][Nb+j-1], 
					e[Nb+i][Nb+j], e[Nb+i][Nb+j+1],
					rhoL, uL, vL, eL, rhoR, uR, vR, eR);
			}
			else if (reconstruction == "MP5_basic") {
				reconstruction_basic_mp5(\
					rho[Nb+i][Nb+j-3], rho[Nb+i][Nb+j-2],
					rho[Nb+i][Nb+j-1], rho[Nb+i][Nb+j],
					rho[Nb+i][Nb+j+1], rho[Nb+i][Nb+j+2],
					u[Nb+i][Nb+j-3], u[Nb+i][Nb+j-2],
					u[Nb+i][Nb+j-1], u[Nb+i][Nb+j],
					u[Nb+i][Nb+j+1], u[Nb+i][Nb+j+2],
					v[Nb+i][Nb+j-3], v[Nb+i][Nb+j-2],
					v[Nb+i][Nb+j-1], v[Nb+i][Nb+j],
					v[Nb+i][Nb+j+1], v[Nb+i][Nb+j+2],
					e[Nb+i][Nb+j-3], e[Nb+i][Nb+j-2],
					e[Nb+i][Nb+j-1], e[Nb+i][Nb+j],
					e[Nb+i][Nb+j+1], e[Nb+i][Nb+j+2],
					rhoL, uL, vL, eL, rhoR, uR, vR, eR);
			}
			else {
				std::cout << "ERROR Fluid2d::rhs > Reconstruction scheme cannnot be specified.";
				exit(EXIT_FAILURE);
			}

			// evaluating rhs using flux scheme
			if (flux_scheme == "Roe_FDS") {
				roe_fds(rhoL, uL, vL, eL, rhoR, uR, vR, eR,
					jxSA, jySA, SA, arr_Fj[i][j]);
			}
			else {
				std::cout << "ERROR Fluid2d::rhs > Flux scheme cannnot be specified.";
				exit(EXIT_FAILURE);
			}
		} // for j
	} // for i


	//--------------------evaluation of RHS---------------------
	for (int i = 0; i < Ni-2*Nb; ++i) {
		for (int j = 0; j < Nj-2*Nb; ++j) {
			for (int k = 0; k < Nf; ++k) {
				arr_Q[i][j][k] = arr_Fi[i][j][k] - arr_Fi[i+1][j][k] \
											 + arr_Fj[i][j][k] - arr_Fj[i][j+1][k]; 
			}
		}
	}


}; // calc_rhs





//================================================================
//=====================Reconstruction Methods=====================
//================================================================

// reconstruct by basic var (rho, u, v, e) with MUSCL-minmod
template<int Ni, int Nj, int Nb>
void Fluid2d<Ni, Nj, Nb>::reconstruction_basic_muscl3(\
	const double& rho0, const double& rho1, 
	const double& rho2, const double& rho3, 
	const double& u0, const double& u1, 
	const double& u2, const double& u3, 
	const double& v0, const double& v1, 
	const double& v2, const double& v3, 
	const double& e0, const double& e1, 
	const double& e2, const double& e3, 
	double& rhoL, double& uL, double& vL, double& eL,
	double& rhoR, double& uR, double& vR, double& eR) const {
	Fnd::MUSCL3(rho0, rho1, rho2, rho3, rhoL, rhoR);
	Fnd::MUSCL3(u0, u1, u2, u3, uL, uR);
	Fnd::MUSCL3(v0, v1, v2, v3, vL, vR);
	Fnd::MUSCL3(e0, e1, e2, e3, eL, eR);
};

// reconstruct by basic var with MP5
template<int Ni, int Nj, int Nb>
void Fluid2d<Ni, Nj, Nb>::reconstruction_basic_mp5(\
	const double& rho0, const double& rho1, 
	const double& rho2, const double& rho3, 
	const double& rho4, const double& rho5, 
	const double& u0, const double& u1, 
	const double& u2, const double& u3,
	const double& u4, const double& u5, 
	const double& v0, const double& v1, 
	const double& v2, const double& v3, 
	const double& v4, const double& v5, 
	const double& e0, const double& e1, 
	const double& e2, const double& e3, 
	const double& e4, const double& e5, 
	double& rhoL, double& uL, double& vL, double& eL,
	double& rhoR, double& uR, double& vR, double& eR) const {
	Fnd::MP5(rho0, rho1, rho2, rho3, rho4, rho5, rhoL, rhoR);
	Fnd::MP5(u0, u1, u2, u3, u4, u5, uL, uR);
	Fnd::MP5(v0, v1, v2, v3, v4, v5, vL, vR);
	Fnd::MP5(e0, e1, e2, e3, e4, e5, eL, eR);
};