#include "fluid2d.h"
//#include<Eigen/Core>
#include <cmath>


// definitions of methods for the derived classes are located
// methods for EoS



// calc eigen values/vectors of flux Jacobian for Euler eq. with ideal gas EoS
template<int Ni, int Nj, int Nb>
void IdealGas2d<Ni, Nj, Nb>::calc_eigen(const double& rho, const double& u, const double& v,
	const double& e, const double& ix, const double& iy,
	double(&LAM)[Nf], double(&R)[Nf][Nf], double(&Rinv)[Nf][Nf]) const {

	// preperation
	const double cs = calc_cs(rho, u, v, e);
	const double p = calc_p(rho, u, v, e);
	const double h = (e + p) / rho;
	const double sqr = std::sqrt(ix * ix + iy * iy);
	const double ixb = ix / sqr;
	const double iyb = iy / sqr;
	const double bigU = ix * u + iy * v;
	const double bigUb = bigU / sqr;
	const double b1 = 0.5 * (u * u + v * v) * (gamma - 1.0) / cs / cs;
	const double b2 = (gamma - 1.0) / cs / cs;

	// calc LAM: diagonal matrix of eigen values
	LAM[0] = bigU - cs * sqr;
	LAM[1] = bigU;
	LAM[2] = bigU + cs * sqr;
	LAM[3] = bigU;

	// calc R: right eigen matrix
	R[0][0] = 1.0;
	R[0][1] = 1.0;
	R[0][2] = 1.0;
	R[0][3] = 0.0;
	R[1][0] = u - ixb * cs;
	R[1][1] = u;
	R[1][2] = u + ixb * cs;
	R[1][3] = -iyb;
	R[2][0] = v - iyb * cs;
	R[2][1] = v;
	R[2][2] = v + iyb * cs;
	R[2][3] = ixb;
	R[3][0] = h - cs * bigUb;
	R[3][1] = 0.5 * (u * u + v * v);
	R[3][2] = h + cs * bigUb;
	R[3][3] = -(iyb * u - ixb * v);

	// calc Rinv: inverse of right eigen matrix
	Rinv[0][0] = 0.5 * (b1 + bigUb / cs);
	Rinv[0][1] = -0.5 * (ixb / cs + b2 * u);
	Rinv[0][2] = -0.5 * (iyb / cs + b2 * v);
	Rinv[0][3] = 0.5 * b2;
	Rinv[1][0] = 1.0 - b1;
	Rinv[1][1] = b2 * u;
	Rinv[1][2] = b2 * v;
	Rinv[1][3] = -b2;
	Rinv[2][0] = 0.5 * (b1 - bigUb / cs);
	Rinv[2][1] = 0.5 * (ixb / cs - b2 * u);
	Rinv[2][2] = 0.5 * (iyb / cs - b2 * v);
	Rinv[2][3] = 0.5 * b2;
	Rinv[3][0] = iyb * u - ixb * v;
	Rinv[3][1] = -iyb;
	Rinv[3][2] = ixb;
	Rinv[3][3] = 0.0;

} // calc_eigen