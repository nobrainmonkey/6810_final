// Filename: main.cpp
// main routine for running a Ising model MCMC chain 
//
// Author: Xihe Han
//
// Changelog:
//	4/3/2023 - initial commit
//	4/7/2023 - added cli menu
//	4/8/2023 - added entropy calculation 
//	4/16/2023 - added microstate graphing mode and gnuplot pipe.
//
//////////////////////////////////////////////////////////////////////////

#include "Observable.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cstdio>

using namespace std;
int main()
{
	// initialize values 
	int thread = 8; // threads to use
	hamiltonian_param_struct hamiltonian_param;
	hamiltonian_param.J = 1.;
	hamiltonian_param.h = 0.01;
	hamiltonian_param_struct* hamiltonian_param_ptr = &hamiltonian_param;
	int row = 4;
	int col = 4;
	int sample = 1000;
	double Tmin = 0.5; //min temperature for evaluation 
	double Tmax = 6; //max temperature for evaluation 
	double dT = 0.1;
	int iteration = row * col * 1000;

	// initial graphing parameters
	double graph_T = 0;
	double graph_time = 0;
	int graph_iteration_skip  = 0;
	int graphing_fps = 10;

	// 	initialize calculated values' vectors 
	std::vector<double> T;
	std::vector<double> Cv;
	std::vector<double> E;
	std::vector<double> Mag;
	std::vector<double> Chi;
	std::vector<double> S;

	// answer to break out of the cli interface
	int answer = 1, answer2 = 1;

	// the cli interface 
	while (answer2 != 0)
	{
		answer = 1;
		while (answer != 0)
		{
			cout << "\nCurrent Calculation Parameters:\n";
			cout << "[1] thread = " << thread << "\t\t";
			cout << "[2] sample = " << sample << "\t" << endl;
			cout << "[3] iterations (auto recommended) = " << iteration << "\t" << endl;
			cout << "[4] J = " << hamiltonian_param.J << "\t";
			cout << "[5] h = " << hamiltonian_param.h << "\t \t" << endl;
			cout << "[6] row = " << row << "\t";
			cout << "[7] col = " << col << endl;
			cout << "[8] Tmin = " << Tmin << "\t";
			cout << "[9] Tmax = " << Tmax << "\t";
			cout << "[10] dT = " << dT << endl;
			cout << "______________________________________"<<endl;
			cout << "Current Graphing Mode Parameters:"<<endl;
			cout << "Parameters cow and rol number are inherited from above:"<<endl;
			cout << "[11] graphing Temperature = " << graph_T <<"\t";
			cout << "[12] graphing time = " << graph_time <<"\t";
			cout << "[13] graphing iteration per frame = " << graph_iteration_skip << "\t" << endl;
			cout << "[14] graphing fps = " << graphing_fps << "\t \t" << endl;
			cout << "[15] start graphing!" << endl;
			cout << "What do you want to change? [0 for none] ";

			cin >> answer;
			cout << endl;

			switch (answer)
			{
				case 0:
					break;
				case 1:
					cout << " enter thread: ";
					cin >> thread;
					break;
				case 2:
					cout << " enter sample size: ";
					cin >> sample;
					break;
				case 3:
					cout << " enter iteration size: ";
					cin >> iteration;
					break;
				case 4:
					cout << " enter J: ";
					cin >> hamiltonian_param.J;
					break;
				case 5:
					cout << " enter h: ";
					cin >> hamiltonian_param.h;
					break;
				case 6:
					cout << " enter row: ";
					cin >> row;
					iteration = row * col * 100;
					break;
				case 7:
					cout << " enter col: ";
					cin >> col;
					iteration = row * col * 100;
					break;
				case 8:
					cout << " enter Tmin: ";
					cin >> Tmin;
					break;
				case 9:
					cout << " enter Tmax: ";
					cin >> Tmax;
					break;

				case 10:
					cout << " enter Temp step ";
					cin >> dT;
					break;
				case 11:
					cout << " enter graphing temp ";
					cin >> graph_T;
					break;
				case 12:
					cout << " enter graphing time ";
					cin >> graph_time;
					break;
				case 13:
					cout << " enter graphing iteration per frame ";
					cin >> graph_iteration_skip;
					break;
				case 14:
					cout << " enter graphing fps ";
					cin >> graphing_fps;
					break;
				case 15:
					cout << " enter to start graphing";
					Microstate *microstate_ptr = new Microstate(row, col, graph_T, hamiltonian_param_ptr);
					microstate_ptr->temperature = graph_T;
					microstate_ptr->graph_evolve(graph_time, graph_iteration_skip, graphing_fps);
					delete microstate_ptr;
					break;
			}
		}

		// clearing calculated values each run
		T.clear();
		E.clear();
		Cv.clear();
		Mag.clear();
		Chi.clear();
		S.clear();

		omp_set_num_threads(thread); // set number of threads to use 

		// output stream for calculated data
		ofstream my_out;
		ostringstream file_name_string;
		file_name_string << "../data/macroscopic_" << setprecision(2) << "J=" << hamiltonian_param.J << "_"
			<< "h=" << hamiltonian_param.h << ".dat";
		string file_name = file_name_string.str();
		my_out.open(file_name.c_str(), ofstream::trunc);

		my_out << "temperature/J"
			<< " " << std::setw(30)
			<< "E/J" << std::setw(37) << " "
			<< "Cv"
			<< " " << std::setw(36)
			<< "m"
			<< " " << std::setw(36)
			<< "chi" << std::endl;

		// stepping through temperatures and gather macroscopic quantities at each temperature.
		for (double temp = Tmin; temp < Tmax; temp += dT)
		{
			std:: cout << "Evaluating Quantities at T = " << temp << std::endl;
			T.push_back(temp);

			double energy = Observable::get_ising_energy(row, col, temp, iteration, sample, &hamiltonian_param);
			E.push_back(energy);

			double cv = Observable::get_ising_heat_capacity(row, col, temp, iteration, sample, &hamiltonian_param);
			Cv.push_back(cv);

			double m = Observable::get_ising_m(row, col, temp, iteration, sample, &hamiltonian_param);
			Mag.push_back(m);

			double chi = Observable::get_ising_chi(row, col, temp, iteration, sample, &hamiltonian_param);
			Chi.push_back(chi);

			my_out << std::scientific << std::setprecision(12)
				<< temp / hamiltonian_param.J << std::setw(20) << " "
				<< energy / hamiltonian_param.J << std::setw(20) << " "
				<< cv << std::setw(20) << " "
				<< m << std::setw(20) << " "
				<< chi << std::setw(20) << " "
				<< std::endl;
		}
		my_out.close();

		// calculate the entropy with the calculated value of Cv
		ostringstream entropy_name_string;
		entropy_name_string << "../data/entropy_" << setprecision(2) << "J=" << hamiltonian_param.J << "_"
			<< "h=" << hamiltonian_param.h << ".dat";
		string entropy_name = entropy_name_string.str();
		my_out.open(entropy_name.c_str(), ofstream::trunc);

		int T_mesh_size = T.size();

		for (int index = 3; index < T_mesh_size; index++) // we start from 3 in order for simpson's rule to work.
		{
			std::vector<double> subT(T.begin(), T.begin() + index+1);
			std::vector<double> subCv(Cv.begin(), Cv.begin() + index+1);
			double entropy = Observable::get_ising_entropy(subT, subCv);
			S.push_back(entropy);
			my_out << subT[index]/hamiltonian_param.J << "    " << entropy << std::endl;
		}
		my_out.close();

		// Creating commands for gnuplot pipe.
		string title = "Macroscopic Quantities";
		string xlabel = "T/k_BJ";
		string ylabel = "J scale";

		stringstream plot_ss;
		stringstream ss;
		ss << "set title '" << title << "'\n";
		ss << "set xlabel '" << xlabel << "'\n";
		ss << "set ylabel '" << ylabel << "'\n";
		ss << "set timestamp " << "\n";
		ss << "set key bottom right " << "\n";
		ss << "plot '" << file_name << "' using 1:2 with lines lw 2 title 'E(T)', \
			'" << file_name << "' using 1:3 with lines lw 2 title 'C_v(T)', \
			'" << entropy_name << "' using 1:2 with lines lw 2 title 'S(T)', \
			'" << file_name << "' using 1:4 with lines lw 2 title 'm(T)', \
			'" << file_name << "' using 1:5 with lines lw 2 title 'chi(T)'\n";
		string plotcmd = ss.str();

		// piping gnuplot command into gnuplot using popen
		FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
		fprintf(gnuplotPipe, "%s", plotcmd.c_str());
		fflush(gnuplotPipe);
		pclose(gnuplotPipe);

		cout << "Again? (no=0, clear=1) ";
		cin >> answer2;
	}
	std::cout << "done!" << std::endl;
	return 0;
}
