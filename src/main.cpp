// include header
#include"fluid2d.h" // declaration of the class IdealGas2d is here
#include<iostream>
#include <string>
#include <time.h>



// include cpp src files
// so that you can compile this project with command: g++ main.cpp
// each file has definitions of methods belong to the class IdealGas2d
#include"fnd.cpp" // reconstruction methods
#include"bc.cpp" // for boundary condition
#include"flux_scheme.cpp" // Euler eq. & flux scheme
#include"marching.cpp" // time-marching scheme & calculation of RHS
#include"io_settings.cpp" // I/O & coordinate-setting
#include"eos.cpp" // EoS




namespace Yahman {
	// estimation of remaining cpu-time
	void estimate_cputime(const int tstep, const int Nout,
		const clock_t cpu_start, float& cpu_time, float& rest_time) {
		// mesure current cpu-time
		clock_t cpu_now = clock();
		// convsrt into sec
		cpu_time = (float)(cpu_now - cpu_start) / CLOCKS_PER_SEC;
		// estimate rest time
		rest_time = cpu_time / (float)tstep * (float)(Nout - tstep);
	}
}



//=====================================================================
//=========================MAIN PROGRAM================================
//=====================================================================
int main()
{

	//--------------------parameter setting------------------------------
	const int Ni = 408; // grids in i-direction
	const int Nj = 408; // grids in j-dirextion
	const int Nb = 4; // number of dummy grids at ends
	const int Nf = 4; // number of variables
	const double Tmax = 3; // time-span to be intrgrated
	const int Nout = 100; // number of file-output
	// for Windows
	//const std::string dir = "..\\output\\";
	//const std::string dir_o = dir + "data_grid408x408\\";
	// for Ubuntu
	const std::string dir = "../output/";
	const std::string dir_o = dir + "data_grid408x408/";
	const std::string f_coordinate = dir + "coordinate_grid408x408.dat";
	const std::string f_settings = dir + "settings.dat";


	//-----------------start cpu-time mesurement--------------------------
	const clock_t cpu_start = clock();



	//--------------------create a fluid object---------------------------
	static IdealGas2d<Ni, Nj, Nb> fluid(dir_o, f_coordinate, f_settings);
	// output settings to external file
	fluid.output_settings(Tmax, Nout);



	// parameters for iteration
	int iter = 0;
	double t = 0.0;
	double dt;
	double Tout = Tmax / (double)Nout;
	float rest_time = 0.0f;
	float cpu_time = 0.0f;



	//-------------------------------iteration----------------------------
	// iteration for output
	for (int tstep = 1; tstep <= Nout; ++tstep) {
		// iteration for marching
		iter = 0;
		while (t < (double)tstep * Tout) {
			// CFL condition
			dt = fluid.calc_cfl(0.7);
			// time marching
			fluid.march_ssprk3(dt, "periodical_in_i", "MP5_basic", "Roe_FDS");
			// update time variables
			t += dt;
			iter += 1;
		}
		// output basic variables
		Yahman::estimate_cputime(tstep, Nout, cpu_start, cpu_time, rest_time);
		fluid.output_basic(tstep, iter, cpu_time, rest_time);
	}


	//---------------------------------END----------------------------------
	std::cout << "Program ended." << std::endl;
} // main