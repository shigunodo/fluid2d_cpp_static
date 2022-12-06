#include "fluid2d.h"
#include <string>
#include<iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
//#include <format>
#include <cstdio>
#include <vector>



// definitions of methods for the class Fluid2d and derived classes are located
// methods for input/output



// format string
namespace Format {
	template <typename ... Args>
	std::string format(const std::string& fmt, Args ... args)
	{
		size_t len = std::snprintf(nullptr, 0, fmt.c_str(), args ...);
		std::vector<char> buf(len + 1);
		std::snprintf(&buf[0], len + 1, fmt.c_str(), args ...);
		return std::string(&buf[0], &buf[0] + len);
	}
}





// input the coordinate from external file
template<int Ni, int Nj, int Nb>
void Fluid2d<Ni, Nj, Nb>::input_coordinate(const std::string& f_name) {

	// open file
	std::ifstream ifs(f_name);
	if (!ifs) {
		std::cout << "ERROR Fluid2d::input_coordinate > Failed to open file." \
			<< std::endl;
		exit(EXIT_FAILURE);
	}

	// input
	for (int i = 0;i < Ni;++i) {
		for (int j = 0;j < Nj;++j) {
			ifs >> x[i][j];
		}
	}
	for (int i = 0;i < Ni;++i) {
		for (int j = 0;j < Nj;++j) {
			ifs >> y[i][j];
		}
	}


	std::cout << "Fluid2d > Coordinate successfully read." << std::endl;
};





// output the coordinate from external file
template<int Ni, int Nj, int Nb>
void Fluid2d<Ni, Nj, Nb>::output_coordinate(const std::string& f_name) const {
	// open file
	std::ofstream ofs(f_name);
	if (!ofs) {
		std::cout << "ERROR Fluid2d::output_coordinate > Failed to open file." \
			<< std::endl;
		exit(EXIT_FAILURE);
	}

	// output
	ofs << std::setprecision(18) << std::scientific;
	for (int i = 0; i < Ni; ++i) {
		for (int j = 0; j < Nj; ++j) {
			ofs << x[i][j] << "\n";
		}
	}
	for (int i = 0; i < Ni; ++i) {
		for (int j = 0; j < Nj; ++j) {
			ofs << y[i][j] << "\n";
		}
	}


	std::cout << "Fluid2d > Coordinate successfully written." << std::endl;
};





// calculate metrices and Jacobian for coordinate, dx for CFL
template<int Ni, int Nj, int Nb>
void Fluid2d<Ni, Nj, Nb>::calc_metrices_dx() {

	// calculate metrices devided by Jacobian
	// evaluated at cell-points (i,j)
	// returning (Ni-2*Nb+2,Nj-2*Nb+2) arrays
	// using 2nd order central difference
	for (int i = 0; i < Ni-2*Nb+2; ++i) {
		for (int j = 0; j < Nj-2*Nb+2; ++j) {
			ixS[i][j] = 0.5 * (y[Nb-1+i][Nb+j] - y[Nb-1+i][Nb-2+j]);
			iyS[i][j] = -0.5 * (x[Nb-1+i][Nb+j] - x[Nb-1+i][Nb-2+j]);
			jxS[i][j] = -0.5 * (y[Nb+i][Nb-1+j] - y[Nb-2+i][Nb-1+j]);
			jyS[i][j] = 0.5 * (x[Nb+i][Nb-1+j] - x[Nb-2+i][Nb-1+j]);
		}
	}

	// calculate inverse of Jacobian
	// evaluated at cell-verteces (i+0.5,j+0.5)
	// returning (Ni-1,Nj-1) array
	for (int i = 0; i < Ni-1; ++i) {
		for (int j = 0; j < Nj-1; ++j) {
			S[i][j] = ((x[1+i][1+j] - x[i][j]) * (y[i][1+j] - y[1+i][j]) \
							- (y[1+i][1+j] - y[i][j]) * (x[i][1+j] - x[1+i][j])) \
							* 0.5;
		}
	}

	// calculate dx for CFL condition
	// evaluated at cell-points (i,j)
	// returning (Ni-2*Nb,Nj-2*Nb) array
	// minimize the adjacent 4 distances
	double dx1;
	double dy1;
	double dd;
	for (int i = 0; i < Ni-2*Nb; ++i) {
		for (int j = 0; j < Nj-2*Nb; ++j) {
			dx1 = x[Nb+1+i][Nb+j] - x[Nb+i][Nb+j];
			dy1 = y[Nb+1+i][Nb+j] - y[Nb+i][Nb+j];
			dx[i][j] = dx1 * dx1 + dy1 * dy1;
			dx1 = x[Nb+i][Nb+1+j] - x[Nb+i][Nb+j];
			dy1 = y[Nb+i][Nb+1+j] - y[Nb+i][Nb+j];
			dd = dx1 * dx1 + dy1 * dy1;
			dx[i][j] = std::min(dx[i][j], dd);
			dx1 = x[Nb-1+i][Nb+j] - x[Nb+i][Nb+j];
			dy1 = y[Nb-1+i][Nb+j] - y[Nb+i][Nb+j];
			dd = dx1 * dx1 + dy1 * dy1;
			dx[i][j] = std::min(dx[i][j], dd);
			dx1 = x[Nb+i][Nb-1+j] - x[Nb+i][Nb+j];
			dy1 = y[Nb+i][Nb-1+j] - y[Nb+i][Nb+j];
			dd = dx1 * dx1 + dy1 * dy1;
			dx[i][j] = std::min(dx[i][j], dd);
			dx[i][j] = std::sqrt(dx[i][j]);
		}
	}
};




// input basic variables from external file
template<int Ni, int Nj, int Nb>
void Fluid2d<Ni, Nj, Nb>::input_basic(const std::string& f_name) {
	// open file
	std::ifstream ifs(f_name);
	if (!ifs) {
		std::cout << "ERROR Fluid2d::input_basic > Failed to open file." \
			<< std::endl;
		exit(EXIT_FAILURE);
	}

	// input
	for (int i = 0;i < Ni;++i) {
		for (int j = 0;j < Nj;++j) {
			ifs >> rho[i][j];
		}
	}
	for (int i = 0;i < Ni;++i) {
		for (int j = 0;j < Nj;++j) {
			ifs >> u[i][j];
		}
	}
	for (int i = 0;i < Ni;++i) {
		for (int j = 0;j < Nj;++j) {
			ifs >> v[i][j];
		}
	}
	for (int i = 0;i < Ni;++i) {
		for (int j = 0;j < Nj;++j) {
			ifs >> e[i][j];
		}
	}

	std::cout << "Fluid2d > Basic variables successfully read." << std::endl;
};





// output basic variables into external file
template<int Ni, int Nj, int Nb>
void Fluid2d<Ni, Nj, Nb>::output_basic(const int tstep, const int iter, 
	float cpu_time, float rest_time) const {
	// open file
	const std::string f_name = Format::format("b%07d.dat", tstep);
	std::ofstream ofs(dir_o + f_name);
	if (!ofs) {
		std::cout << "ERROR Fluid2d::output_basic > Failed to open file." \
			<< std::endl;
		exit(EXIT_FAILURE);
	}

	// write basic variables into file
	ofs << std::setprecision(18) << std::scientific;
	for (int i = 0;i < Ni;++i) {
		for (int j = 0;j < Nj;++j) {
			ofs << rho[i][j] << "\n";
		}
	}
	for (int i = 0;i < Ni;++i) {
		for (int j = 0;j < Nj;++j) {
			ofs << u[i][j] << "\n";
		}
	}
	for (int i = 0;i < Ni;++i) {
		for (int j = 0;j < Nj;++j) {
			ofs << v[i][j] << "\n";
		}
	}
	for (int i = 0;i < Ni;++i) {
		for (int j = 0;j < Nj;++j) {
			ofs << e[i][j] << "\n";
		}
	}

	// elapsed time convertion
	constexpr float sec_per_min = 60.0f;
	int hh = (int)(cpu_time / (sec_per_min * sec_per_min));
	cpu_time = cpu_time - (float)hh * sec_per_min * sec_per_min;
	int mm = (int)(cpu_time / sec_per_min);
	cpu_time = cpu_time - (float)mm * sec_per_min;
	int ss = (int)cpu_time;

	// rest time convertion
	int hh_r = (int)(rest_time / (sec_per_min * sec_per_min));
	rest_time = rest_time - (float)hh_r * sec_per_min * sec_per_min;
	int mm_r = (int)(rest_time / sec_per_min);
	rest_time = rest_time - (float)mm_r * sec_per_min;
	int ss_r = (int)rest_time;

	// output status to standard
	std::cout << "elapsed: " << std::setfill(' ') << std::setw(3) << hh << " h " \
		<< std::setfill('0') << std::setw(2) << mm << " m " \
		<< std::setfill('0') << std::setw(2) << ss << " s | " \
		<< "tstep = " << std::setfill(' ') << std::setw(5) << tstep << " | " \
		<< "t = " << std::setprecision(4) << std::setw(10) << t << " | " \
		<< "iter = " << std::setw(7) << iter << " | " \
		<< "rest: " << std::setfill(' ') << std::setw(3) << hh_r << " h " \
		<< std::setfill('0') << std::setw(2) << mm_r << " m " \
		<< std::setfill('0') << std::setw(2) << ss_r << " s " \
		<< std::setfill(' ') << std::endl;

	// output status into settings-file
	std::ofstream ofs_s;
	ofs_s.open(f_settings, std::ios::app);
	ofs_s << "elapsed: " << std::setfill(' ') << std::setw(3) << hh << " h " \
		<< std::setfill('0') << std::setw(2) << mm << " m " \
		<< std::setfill('0') << std::setw(2) << ss << " s | " \
		<< "tstep = " << std::setfill(' ') << std::setw(5) << tstep << " | " \
		<< "t = " << std::setprecision(4) << std::setw(10) << t << " | " \
		<< "iter = " << std::setw(7) << iter << " | " \
		<< "rest: " << std::setfill(' ') << std::setw(3) << hh_r << " h " \
		<< std::setfill('0') << std::setw(2) << mm_r << " m " \
		<< std::setfill('0') << std::setw(2) << ss_r << " s " \
		<< std::setfill(' ') << std::endl;
	ofs_s.close();

};





// output settings, for class IdealGas2d
template<int Ni, int Nj, int Nb>
void IdealGas2d<Ni, Nj, Nb>::output_settings(const double Tmax, const int Nout) {
	// open file
	std::ofstream ofs(this->f_settings);
	if (!ofs) {
		std::cout << "ERROR IdealGas2d::output_settings > Failed to open file." \
			<< std::endl;
		exit(EXIT_FAILURE);
	}
	ofs << "Solving Euler equation system with ideal gas eos. \n" \
		<< "directory: " << this->dir_o << "\n" \
		<< "gamma = " << std::setprecision(4) << std::setw(10) << gamma << "\n" \
		<< "Ni    = " << std::setfill(' ') << std::setw(10) << Ni << "\n" \
		<< "Nj    = " << std::setfill(' ') << std::setw(10) << Nj << "\n" \
		<< "Nb    = " << std::setfill(' ') << std::setw(10) << Nb << "\n" \
		<< "Tmax  = " << std::setprecision(4) << std::setw(10) << Tmax << "\n" \
		<< "Nout  = " << std::setfill(' ') << std::setw(10) << Nout << "\n" \
		<< "\n";
		
	std::cout << "IdealGas2d > Calculation settings file outputted." \
		<< std::endl;
};