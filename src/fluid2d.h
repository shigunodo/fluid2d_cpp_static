#pragma once
#include <string>
#include <cmath>
#include <time.h>
#include<iostream>

// Contents:
// + declaration of virtual class Fluid2d, which diesn't have the definition of EoS
// + declaration of class IdealGas2d derived from Fluid2d, which has ideal EoS



//==================================================================
//--------------------Base Virtual Class Fluid2d--------------------
// without EoS definition
//==================================================================
template<int Ni, int Nj, int Nb>// templete patameter
class Fluid2d
{


	//--------------------------------------------------------------
	//========================Protected Members=====================
	//--------------------------------------------------------------
protected:



	//-----------------------Global Variables-----------------------

	// parameters

	static const int Nf = 4; // size of flux as vector

	// coordinate

	double t = 0.0; // calculation time
	int tstep = 0; // time steps

	double x[Ni][Nj]; // x-coordinate list of grids: (Ni,Nj) array, evaluated at (i,j)
	double y[Ni][Nj]; // y-coordinate list of grids: (Ni,Nj) array, evaluated at (i,j)

	// basic variables of fluid

	double rho[Ni][Nj]; // density: (Ni,Nj) array, evaluated at (i,j)
	double u[Ni][Nj]; // x-component of velocity: (Ni,Nj) array, evaluated at (i,j)
	double v[Ni][Nj]; // y-component of velocity: (Ni,Nj) array, evaluated at (i,j)
	double e[Ni][Nj]; // total energy density per volume: (Ni,Nj) array, evaluated at (i,j)
	
	// temporal variables needed for methods

	double arr_Q0[Ni-2*Nb][Nj-2*Nb][Nf]; // for march_ssprk3
	double arr_Q1[Ni-2*Nb][Nj-2*Nb][Nf]; // for march_ssprk3
	double arr_Q2[Ni-2*Nb][Nj-2*Nb][Nf]; // for march_ssprk3

	double arr_Fi[Ni-2*Nb+1][Nj-2*Nb][Nf]; // for rhs
	double arr_Fj[Ni-2*Nb][Nj-2*Nb+1][Nf]; // for rhs

	// directory for data-output
	const std::string dir_o;
	
	// file name for coordinate
	const std::string f_coordinate;

	// file name for settings-file
	const std::string f_settings;

	// metrices and Jacobian
	
	// inverse of Jacobian: (Ni-1,Nj-1) array, evaluated at (i+0.5,j+0.5)
	double S[Ni-1][Nj-1];
	// metrix ix times S: (Ni-2*Nb+2,Nj-2*Nb+2) array, evaluated at (i,j)
	double ixS[Ni-2*Nb+2][Nj-2*Nb+2];
	// metrix iy times S: (Ni-2*Nb+2,Nj-2*Nb+2) array, evaluated at (i,j)
	double iyS[Ni-2*Nb+2][Nj-2*Nb+2];
	// metrix jx times S: (Ni-2*Nb+2,Nj-2*Nb+2) array, evaluated at (i,j)
	double jxS[Ni-2*Nb+2][Nj-2*Nb+2];
	// metrix jy times S: (Ni-2*Nb+2,Nj-2*Nb+2) array, evaluated at (i,j)
	double jyS[Ni-2*Nb+2][Nj-2*Nb+2];

	// dx for CFL condition: (Ni-2*Nb,Nj-2*Nb) array, evaluated at (i,j)
	double dx[Ni-2*Nb][Nj-2*Nb];




	//-------Virtual Functions Need to be Defined in Derived Classes-------

	// about equation of state

	// calc pressure from density, velocity, total energy per volume
	virtual double calc_p(double rho, double u, double v, double e) const = 0;

	// calc temperature from density, velocity, total energy per volume
	virtual double calc_T(double rho, double u, double v, double e) const = 0;

	// calc sound-speed from density, velocity, total energy per volume
	virtual double calc_cs(double rho, double u, double v, double e) const = 0;

	// calc total energy per volume
	// from density, velocity, total specific enthalpy
	virtual double calc_e(double rho, double u, double v, double h) const = 0;

	// calc density and total energy per volume
	// from pressure, temperature, velocity
	virtual void calc_rho_e(const double& p, const double& T, const double& u,
		const double& v, double& rho, double& e) const = 0;

	// calc total energy per volume from density, velocity, pressure
	virtual double calc_e_wp(double rho, double u, double v, double p) const = 0;

	// eigen values/vectors of flux Jacobian
	virtual void calc_eigen(const double& rho, const double& u, const double& v,
		const double& e, const double& ix, const double& iy,
		double(&LAM)[Nf], double(&R)[Nf][Nf], double(&Rinv)[Nf][Nf]) const = 0;





	//-------------------Convenient Functions----------------------

	// calc velocity norm from components
	double calc_v(const double u, const double v) const
	{
		return std::sqrt(u * u + v * v);
	};
	



	//---------Methods for Flux/conserved Construction----------
	// definitions are located in "flux_scheme.cpp"

	// construct conserved var from basic var
	void calc_conserved(const double& rho, const double& u,
		const double& v, const double& e, const double& S,
		double(&Q)[Nf]) const;

	// construct convection fux from basic var
	void calc_flux_conv(const double& rho, const double& u,
		const double& v, const double& e, const double& ixS,
		const double& iyS, double(&Fc)[Nf]) const;

	// construct basic var from conserved var
	void calc_basic(const double(&Q)[Nf], const double& S,
		double& rho, double& u, double& v, double& e) const;




	//------------------------Flux Schemes-------------------------
	// definitions are located in "flux_scheme.cpp"

	// calculate Roe average
	void roe_average(const double& rhoL, const double& uL, const double& vL,
		const double& eL, const double& rhoR, const double& uR,
		const double& vR, const double& eR, double& rhoA,
		double& uA, double& vA, double& eA) const;
	//
	void roe_fds(const double& rhoL, const double& uL, const double& vL,
		const double& eL, const double& rhoR, const double& uR,
		const double& vR, const double& eR, const double& ixS,
		const double& iyS, const double& S, double(&Fc)[Nf]) const;





	//----------------Evaluating RHS of Equation-------------------
	void rhs(const double(&rho)[Ni][Nj], const double(&u)[Ni][Nj],
		const double(&v)[Ni][Nj], const double(&e)[Ni][Nj],
		double(&arr_Q)[Ni-2*Nb][Nj-2*Nb][Nf], 
		const std::string& reconstruction = "MP5_basic",
		const std::string& flux_scheme = "Roe_FDS");

	// reconstruction schemes
	void reconstruction_basic_muscl3(\
		const double& rho0, const double& rho1, 
		const double& rho2, const double& rho3, 
		const double& u0, const double& u1, 
		const double& u2, const double& u3, 
		const double& v0, const double& v1, 
		const double& v2, const double& v3, 
		const double& e0, const double& e1, 
		const double& e2, const double& e3, 
		double& rhoL, double& uL, double& vL, double& eL,
		double& rhoR, double& uR, double& vR, double& eR) const;
	void reconstruction_basic_mp5(\
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
		double& rhoR, double& uR, double& vR, double& eR) const;




	//---------------------Boundary Conditions---------------------
	// definitions are located in "bc.cpp"

	// reflecting Boundary Condition of selected type
	void reflect_bc(const std::string& bc_type);

	// periodical in i-direction, Dirichlet on j-boundaries
	void bc_periodical_in_i();
	




	//-------------------------------------------------------------
	//========================Public Members=======================
	//-------------------------------------------------------------
public:


	//-----------------------Constructor---------------------------
	explicit Fluid2d( \
		const std::string& dir_o, 
		const std::string& f_coordinate,
		const std::string& f_settings)
		: dir_o(dir_o), f_coordinate(f_coordinate), f_settings(f_settings)
	{	
		// initialize: input coordinate and initial condition,
		// calculate metrices, Jacobian and dx for CFL
		initialize();
	};




	//------------------time marching schemes----------------------
	// definitions are located in "marching.cpp"

	// 3rd order SSP Runge-Kutta method
	void march_ssprk3(const double& dt, 
		const std::string& bc_type,
		const std::string& reconstruction = "MP5_basic",
		const std::string& flux_scheme = "Roe_FDS");




	//----------------------CFL condition--------------------------
	// definitions are located in "marching.cpp"
	double calc_cfl(double CFL_coeff) const;




	//----------Methods for i/o and Initial Settings---------------
	// definitions are located in "io_settings.cpp"

	// input the coordinate from external file
	void input_coordinate(const std::string& f_name);

	// output the coordinate from external file
	void output_coordinate(const std::string& f_name) const;

	// calculate metrices and Jacobian for coordinate, dx for CFL
	void calc_metrices_dx();

	// input basic variables from external file
	void input_basic(const std::string& f_name);

	// output basic variables into external file
	void output_basic(int tstep, int iter, 
		float cpu_time, float rest_time = 0.0f) const;

	// initialize fluid: setting coordinate and inputting initial condition
	void initialize() {
		t = 0.0;
		tstep = 0;
		input_coordinate(f_coordinate);
		input_basic(dir_o + "b0000000.dat");
		calc_metrices_dx();
	};




	//--------------------Other Methods----------------------------
	double get_time() const {
		return t;
	};

	int get_tstep() const {
		return tstep;
	};


}; // class Fluid2d






//========================================================================
//-----------------Derived Class with Ideal Gas EOS-----------------------
//========================================================================
template<int Ni, int Nj, int Nb>
class IdealGas2d : public Fluid2d<Ni, Nj, Nb>
{
protected:

	static const int Nf = 4; // size of flux as vector

	// specific heat
	const double gamma;




	//---------------calculation methods for ideal gas EoS----------------

	// get the specific heat
	double get_gamma() const { return gamma; };

	// calc pressure from density, velocity, total energy per volume
	double calc_p(const double rho, const double u,
		const double v, const double e) const override {
		return (gamma - 1.0) * (e - 0.5 * rho * (u * u + v * v));
	};

	// calc temperature from density, velocity, total energy per volume
	double calc_T(const double rho, const double u,
		const double v, const double e) const override {
		return e / rho - 0.5 * (u * u + v * v);
	};

	// calc sound-speed from density, velocity, total energy per volume
	double calc_cs(const double rho, const double u,
		const double v, const double e) const override {
		return std::sqrt(gamma * (gamma - 1.0) * (e / rho - 0.5 * (u * u + v * v)));
	};

	// calc total energy per volume
	// from density, velocity, total specific enthalpy
	double calc_e(const double rho, const double u,
		const double v, const double h) const override {
		return rho * (h + (gamma - 1.0) * 0.5 * (u * u + v * v)) / gamma;
	};

	// calc density and total energy per volume
	// from pressure, temperature, velocity
	void calc_rho_e(const double& p, const double& T, const double& u,
		const double& v, double& rho, double& e) const override {
		rho = p / T / (gamma - 1.0);
		e = rho * (T + 0.5 * (u * u + v * v));
	};

	// calc total energy per volume from density, velocity, pressure
	double calc_e_wp(const double rho, const double u,
		const double v, const double p) const override {
		return p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v);
	};

	// eigen values/vectors of flux Jacobian
	// definition is located in "eos.cpp"
	void calc_eigen(const double& rho, const double& u, const double& v,
		const double& e, const double& ix, const double& iy,
		double(&LAM)[Nf], double(&R)[Nf][Nf], double(&Rinv)[Nf][Nf]) const override;



public:



	//----------------------------Constructor--------------------------
	explicit IdealGas2d( \
		const std::string& dir_o,
		const std::string& f_coordinate,
		const std::string& f_settings,
		const double gamma = 1.4)
		: Fluid2d<Ni, Nj, Nb>(dir_o, f_coordinate, f_settings), 
		gamma(gamma) {};




	//------------------------------I/O--------------------------------

	// output settings
	// definition is located in "io_settings.cpp"
	void output_settings(double Tmax, int Nout);
};

