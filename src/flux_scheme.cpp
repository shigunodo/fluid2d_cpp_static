#include "fluid2d.h"
#include <cmath>


// definitions of methods for the class Fluid2d are located
// methods for flux schemes




//---------Methods for Flux/conserved Construction----------




// construct conserved var from basic var
template<int Ni, int Nj, int Nb>
void Fluid2d<Ni, Nj, Nb>::calc_conserved(const double& rho, const double& u,
	const double& v, const double& e, const double& S, 
	double(&Q)[Nf]) const {
	Q[0] = rho * S;
	Q[1] = rho * u * S;
	Q[2] = rho * v * S;
	Q[3] = e * S;
};



// construct convection flux from basic var
template<int Ni, int Nj, int Nb>
void Fluid2d<Ni, Nj, Nb>::calc_flux_conv(const double& rho, const double& u,
	const double& v, const double& e, const double& ixS,
	const double& iyS, double(&Fc)[Nf]) const {
	// preperation
	const double bigU = ixS * u + iyS * v;
	const double p = calc_p(rho, u, v, e);
	Fc[0] = rho * bigU;
	Fc[1] = rho * u * bigU + ixS * p;
	Fc[2] = rho * v * bigU + iyS * p;
	Fc[3] = (e + p) * bigU;
};



// construct basic var from conserved var
template<int Ni, int Nj, int Nb>
void Fluid2d<Ni, Nj, Nb>::calc_basic(const double(&Q)[Nf], const double& S,
	double& rho, double& u, double& v, double& e) const {
	rho = Q[0] / S;
	u = Q[1] / rho / S;
	v = Q[2] / rho / S;
	e = Q[3] / S;
};





//------------------------Flux Schemes-------------------------




// calculate Roe average from L/R value
template<int Ni, int Nj, int Nb>
void Fluid2d<Ni, Nj, Nb>::roe_average(const double& rhoL, const double& uL, const double& vL,
	const double& eL, const double& rhoR, const double& uR,
	const double& vR, const double& eR, double& rhoA,
	double& uA, double& vA, double& eA) const {

	// preperation
	const double pL = calc_p(rhoL, uL, vL, eL);
	const double pR = calc_p(rhoR, uR, vR, eR);
	const double hL = (eL + pL) / rhoL;
	const double hR = (eR + pR) / rhoR;
	const double sqrL = std::sqrt(rhoL);
	const double sqrR = std::sqrt(rhoR);

	// averaging
	rhoA = sqrL * sqrR;
	uA = (sqrL * uL + sqrR * uR) / (sqrL + sqrR);
	vA = (sqrL * vL + sqrR * vR) / (sqrL + sqrR);
	const double hA = (sqrL * hL + sqrR * hR) / (sqrL + sqrR);
	eA = calc_e(rhoA, uA, vA, hA);
};



// classical Roe-type FDS scheme (Roe, 1981)
template<int Ni, int Nj, int Nb>
void Fluid2d<Ni, Nj, Nb>::roe_fds(const double& rhoL, const double& uL, const double& vL,
	const double& eL, const double& rhoR, const double& uR,
	const double& vR, const double& eR, const double& ixS,
	const double& iyS, const double& S, double(&Fc)[Nf]) const {

	//preperation
	const double ix = ixS / S;
	const double iy = iyS / S;

	// construction of some matrices
	double rhoA, uA, vA, eA;
	roe_average(rhoL, uL, vL, eL, rhoR, uR, vR, eR, rhoA, uA, vA, eA);
	double QL[Nf], QR[Nf];
	calc_conserved(rhoL, uL, vL, eL, S, QL);
	calc_conserved(rhoR, uR, vR, eR, S, QR);
	double FL[Nf], FR[Nf];
	calc_flux_conv(rhoL, uL, vL, eL, ixS, iyS, FL);
	calc_flux_conv(rhoR, uR, vR, eR, ixS, iyS, FR);
	double LAM[Nf], R[Nf][Nf], Rinv[Nf][Nf];
	calc_eigen(rhoA, uA, vA, eA, ix, iy, LAM, R, Rinv);

	// calc numerical flux at cell-boundaries by FDS formulation
	for (int i = 0; i < Nf; ++i) {
		Fc[i] = FR[i] + FL[i];
		for (int j = 0; j < Nf; ++j) {
			for (int k = 0; k < Nf; ++k) {
				Fc[i] -= R[i][j] * std::abs(LAM[j]) * Rinv[j][k] \
																						* (QR[k] - QL[k]);
			}
		}
		Fc[i] *= 0.5;
	}


}; // roe_fds