#pragma once
#include <algorithm>
#include <cmath>


// declaration of class Fnd, having methods for general purposes



//-----------------methods for reconstruction schemes----------------
class Fnd
{
private:
	// ingredients of each scheme

	// returning sign of arg
	static double sign(const double x)
	{
		return (double)(x > 0.0) - (double)(x < 0.0);
	};



	// minmod limitter for MUSCL and MP5
	static double minmod(const double x, const double y)
	{
		double sgn = sign(x);
		return sgn * std::max(std::min(std::abs(x), sgn * y), 0.0);
	};



	// minmod of 4 doubles for MP5
	static double minmod4(const double x1, const double x2, 
		const double x3, const double x4)
	{
		return 0.5 * (sign(x1) + sign(x2)) * std::abs(\
			0.5 * (sign(x1) + sign(x3)) \
			* 0.5 * (sign(x1) + sign(x4))) \
			* std::min({ std::abs(x1),std::abs(x2),std::abs(x3),std::abs(x4) });
	};



	// median for MP5
	static double median(const double x, const double y, 
		const double z)
	{
		return x + minmod(y - x, z - x);
	};



public:



	//-------------------Reconstruction Scheme------------------
	// definitions are located in "fnd.cpp"


	// MUSCL: k=1/3 (van Leer, 1977) w/ minmod limitter
	// evaluating L/R at cell-boundary between q and qp
	static void MUSCL3(const double& qm, const double& q,
		const double& qp, const double& q2p,
		double& qL, double& qR);


	// MP5: alp=2 (Suresh & Huynh, 1997)
	// evaluating L/R at cell-boundary between q and qp
	static void MP5_sub(const double& q2m, const double& qm,
		const double& q, const double& qp,
		const double& q2p, double& qL);


	static void MP5(const double& q2m, const double& qm,
		const double& q, const double& qp,
		const double& q2p, const double q3p,
		double& qL, double& qR)
	{
		MP5_sub(q2m, qm, q, qp, q2p, qL);
		MP5_sub(q3p, q2p, qp, q, qm, qR);
	};


}; // class Fnd

