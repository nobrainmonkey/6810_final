#include "Observable.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <string>
#include <sstream>

using namespace std;
int main()
{
	int thread = 8;
	hamiltonian_param_struct hamiltonian_param;
	hamiltonian_param.J = 1.;
	hamiltonian_param.h = 0.;
	int row = 20;
	int col = 20;
	int sample = 200;
	double Tmin = 0.5;
	double Tmax = 6;
	double dT = 0.1;
	std::vector<double> T;
	std::vector<double> Cv;
	std::vector<double> E;
	std::vector<double> Mag;
	std::vector<double> Chi;
	std::vector<double> S;

	int answer = 1, answer2 = 1;
	while (answer2 != 0)
	{
		answer = 1;
		while (answer != 0)
		{
			cout << "\nCurrent parameters:\n";
			cout << "[1] thread = " << thread << "\t\t";
			cout << "[2] sample = " << sample << "\t" << endl;
			cout << "[3] J = " << hamiltonian_param.J << "\t";
			cout << "[4] h = " << hamiltonian_param.h << "\t" << endl;
			cout << "[5] row = " << row << "\t";
			cout << "[6] col = " << col << endl;
			cout << "[7] Tmin = " << Tmin << "\t";
			cout << "[8] Tmax = " << Tmax << "\t";
			cout << "[9] dT = " << dT << endl;
			cout << "\nWhat do you want to change? [0 for none] ";

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
				cout << " enter J: ";
				cin >> hamiltonian_param.J;
				break;
			case 4:
				cout << " enter h: ";
				cin >> hamiltonian_param.h;
				break;
			case 5:
				cout << " enter row: ";
				cin >> row;
				break;
			case 6:
				cout << " enter col: ";
				cin >> col;
				break;
			case 7:
				cout << " enter Tmin: ";
				cin >> Tmin;
				break;
			case 8:
				cout << " enter Tmax: ";
				cin >> Tmax;
				break;

			case 9:
				cout << " enter Temp step ";
				cin >> dT;
				break;
			}
		}
		T.clear();
		E.clear();
		Cv.clear();
		Mag.clear();
		Chi.clear();
		S.clear();
		
		omp_set_num_threads(thread);
		int iteration = row * col * 1000;
		ofstream my_out;
		ostringstream file_name_string;
		file_name_string << "../data/macroscopic_" << setprecision(2) << "J=" << hamiltonian_param.J << "_"
						 << "h=" << hamiltonian_param.h << ".dat";
		string file_name = file_name_string.str();
		my_out.open(file_name.c_str(), ofstream::trunc);

		my_out << "tempreature/J"
			   << " " << std::setw(30)
			   << "E/J" << std::setw(37) << " "
			   << "Cv"
			   << " " << std::setw(36)
			   << "m"
			   << " " << std::setw(36)
			   << "chi" << std::endl;

		for (double temp = Tmin; temp < Tmax; temp += dT)
		{
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

		ostringstream entropy_name_string;
		entropy_name_string << "../data/entropy_" << setprecision(2) << "J=" << hamiltonian_param.J << "_"
							<< "h=" << hamiltonian_param.h << ".dat";
		string entropy_name = entropy_name_string.str();
		my_out.open(entropy_name.c_str(), ofstream::trunc);

		int T_mesh_size = T.size();

		for (int index = 3; index < T_mesh_size; index++)
		{
			std::vector<double> subT(T.begin(), T.begin() + index);
			std::vector<double> subCv(Cv.begin(), Cv.begin() + index);
			double entropy = Observable::get_ising_entropy(subT, subCv);
			S.push_back(entropy);
			my_out << subT[index]/hamiltonian_param.J << "    " << entropy << std::endl;
		}
		my_out.close();

		cout << "Again? (no=0, clear=1) ";
		cin >> answer2;
	}
	std::cout << "done!" << std::endl;
	return 0;
}
