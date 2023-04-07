// Fila name: Observable.cpp
//
// Source file for calculating observables from a microstate
//
// Programmer: Xihe Han
//
// Version history:
//	04/06/2023: original version
//
//**************************************************************

// include files
#include "Hamiltonian.h"
#include "Microstate.h"
#include "Observable.h"
#include <cmath>
#include <iostream>
#include <omp.h>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

double spline_helper(double x, void* params)
{

}

// return the ising energy of the single spin by doing the folowing:
// 1. initialize a microstate, evolve the microstate for `iteration` times so it reaches thermal equalibrium
// 2. keep evolving the microstate for `mc_steps` amount of time, record the energy for each evolution, and average the energy
// 3. return to step one for `sample_size` amount of time, record the averaged energy each time, and finally average all the samples
// 4. return the total sample energy per spin in the unit of J.
double Observable::get_ising_energy(int row, int col, double T, int iteration, int sample_size, hamiltonian_param_struct *hamiltonian_param_ptr){
	double total_sample_energy = 0.;
#pragma omp parallel for
	for(int i=0; i<sample_size; i++)
	{
		Microstate* microstate_ptr = new Microstate(row, col, T, hamiltonian_param_ptr);

		// evolve the microstate to reach equalibruim
		microstate_ptr->evolve_microstate(iteration);

		// create a list to store energy after some iterations 
		// we choose a list because it is dynamically sized
		// and thread safe in case we ant to use omp.

		// after reaching equalibrium, we want to https://discordapp.com/channels/@me/1028655396090544259 further
		// evolve it 1000 time and record energy for each time 
		// and take the final average.
		int mc_steps = 1000;    // amount of step to take average
		// get the total amount of energy over the num_average iterations of the microstate.
		double total_energy = 0;
		for(int i=0; i < mc_steps; i++){
			microstate_ptr->evolve_microstate(1);
			total_energy += (Hamiltonian::hamiltonian_periodic_ising_microstate(microstate_ptr->get_microstate_matrix_ptr(), hamiltonian_param_ptr));
		}
		double microstate_total_energy = total_energy / double(mc_steps);
		double microstate_spin_energy = microstate_total_energy / (hamiltonian_param_ptr->J * double((row * col)));
		total_sample_energy += microstate_spin_energy;
		delete microstate_ptr;
	}
	return total_sample_energy/double(sample_size);
}

// return the ising heat capacity of the single spin by doing the folowing:
// 1. initialize a microstate, evolve the microstate for `iteration` times so it reaches thermal equalibrium
// 2. keep evolving the microstate for `mc_steps` amount of time, record the energy for each evolution, and use Cv = (<E^2> - <E>^2)/T to calcualte C_v
// 3. return to step one for `sample_size` amount of time, record the averaged energy each time, calculate Cv, and finally average all the samples
// 4. return the total sample specific heat per spin in the unit of J.
double Observable::get_ising_heat_capacity(int row, int col, double T, int iteration, int sample_size, hamiltonian_param_struct* hamiltonian_param_ptr)
{
	double total_sample_heat_capacity = 0;
#pragma omp parallel for
	for (int i = 0; i < sample_size; i++)
	{
		Microstate* microstate_ptr = new Microstate(row, col, T, hamiltonian_param_ptr);
		microstate_ptr->evolve_microstate(iteration);    //bring system into thermal equalibrium
		double E = 0.;
		double Esq = 0.;
		int mc_steps = 1000;    // steps after thermal equalibrium
		for(int i = 0; i < mc_steps; i++){
			microstate_ptr->evolve_microstate(1);
			double energy = Hamiltonian::hamiltonian_periodic_ising_microstate(microstate_ptr->get_microstate_matrix_ptr(), hamiltonian_param_ptr);
			E += energy;
			Esq += energy * energy;
		}
		double E_avg = E / mc_steps;
		double Esq_avg = Esq / mc_steps;
		double specific_heat = (Esq_avg - E_avg * E_avg) / (T * T * double(row * col)* hamiltonian_param_ptr->J);
		total_sample_heat_capacity += specific_heat;
		delete microstate_ptr;
	}
	return total_sample_heat_capacity / (double)sample_size;

}

// return the ising entropy per particle per J of the single spin by using the equation
// S = U/T + ln(Cv). Note the user should use the respective method to calculate E and Cv.
double Observable::get_ising_entropy(int row, int col, int iteration, int sample_size, hamiltonian_param_struct* hamiltonian_param_ptr, double lower_T, double upper_T, double dT){
	double entropy = 0.;
	double cv_prev = get_ising_heat_capacity(row, col, lower_T, iteration, sample_size, hamiltonian_param_ptr);
	double cv_curr = get_ising_heat_capacity(row, col, lower_T+dT, iteration, sample_size, hamiltonian_param_ptr);
	double cv_next;
	for (double T = lower_T + 2.* dT ; T < upper_T; T += 2. * dT) {
		cv_next = get_ising_heat_capacity(row, col, T, iteration, sample_size, hamiltonian_param_ptr);
		entropy += (dT/3.0) * (cv_prev/(T-2*dT) + 4*cv_curr/(T-dT) + cv_next/T);
		cv_prev = cv_curr;
		cv_curr = cv_next;
}

    return entropy;

}

double Observable::get_ising_m(int row, int col, double T, int iteration, int sample_size, hamiltonian_param_struct* hamiltonian_param_ptr)
{
	double total_sample_m = 0.;
#pragma omp parallel for
	for(int i=0; i<sample_size; i++)
	{
		Microstate* microstate_ptr = new Microstate(row, col, T, hamiltonian_param_ptr);

		// evolve the microstate to reach equalibruim
		microstate_ptr->evolve_microstate(iteration);

		// create a list to store energy after some iterations 
		// we choose a list because it is dynamically sized
		// and thread safe in case we ant to use omp.

		// after reaching equalibrium, we want to https://discordapp.com/channels/@me/1028655396090544259 further
		// evolve it 1000 time and record energy for each time 
		// and take the final average.
		int mc_steps = 1000;    // amount of step to take average
		// get the total amount of energy over the num_average iterations of the microstate.
		double total_m = 0;
		for(int i=0; i < mc_steps; i++){
			microstate_ptr->evolve_microstate(1);
			total_m += std::fabs((microstate_ptr->get_microstate_matrix_ptr())->sum());
		}
		double average_m = total_m / mc_steps;
		total_sample_m += average_m / double(row * col);
		delete microstate_ptr;
	}
	return total_sample_m / double(sample_size);
}

double Observable::get_ising_chi(int row, int col, double T, int iteration, int sample_size, hamiltonian_param_struct *hamiltonian_param_ptr)
{
	double total_sample_chi = 0;
#pragma omp parallel for
	for (int i = 0; i < sample_size; i++)
	{
		Microstate* microstate_ptr = new Microstate(row, col, T, hamiltonian_param_ptr);
		microstate_ptr->evolve_microstate(iteration);    //bring system into thermal equalibrium
		double m = 0.;
		double msq = 0.;
		int mc_steps = 1000;    // steps after thermal equalibrium
		for(int i = 0; i < mc_steps; i++){
			microstate_ptr->evolve_microstate(1);
			double magnetization = microstate_ptr->get_microstate_matrix_ptr()->mean();
			m += magnetization;
			msq += magnetization * magnetization;
		}
		double m_avg = m / mc_steps;
		double msq_avg = msq / mc_steps;
		double chi = (msq_avg - m_avg * m_avg) / (T * double(row * col)* hamiltonian_param_ptr->J);
		total_sample_chi += chi;
		delete microstate_ptr;
	}
	return total_sample_chi / (double)sample_size;
}



int main()
{
	omp_set_num_threads(16);
	hamiltonian_param_struct hamiltonian_param;
	hamiltonian_param.J = 2.;
	hamiltonian_param.h = 0.1;
	int row = 4;
	int col = 4;
	int iteration = 1024;
	int sample = iteration;
	double dT = 0.01;

	std::ofstream my_out;
	my_out.open("output.dat");
	my_out << "tempreature" << " " << std::setw(28) 
		   << "E"  << std::setw(37) << " "
		   << "Cv" << " " << std::setw(36)
//		   << "S" << " " << std::setw(36)
		   << "m" << " " << std::setw(36)
		   << "chi" << std::endl;
	for (double temp = 1.5; temp < 15.; temp += 0.1)
	{
		double energy = Observable::get_ising_energy(row, col, temp, iteration, sample, &hamiltonian_param);
		double cv = Observable::get_ising_heat_capacity(row, col, temp, iteration, sample, &hamiltonian_param);
//		double entropy = Observable::get_ising_entropy(row, col, iteration, sample, &hamiltonian_param, 2. * dT, temp, dT); 
		double m = Observable::get_ising_m(row, col, temp, iteration, sample, &hamiltonian_param);
		double chi = Observable::get_ising_chi(row, col, temp, iteration, sample, &hamiltonian_param);
		my_out <<std::scientific << std::setprecision(12)
			  << temp/hamiltonian_param.J << std::setw(20) <<" " 
			  << energy <<  std::setw(20) << " "
			  << cv  << std::setw(20) << " "
//			  << entropy <<std::setw(20) << " "
			  << m <<std::setw(20) << " "
			  << chi <<std::setw(20) << " "
			  <<std::endl; 
	}

	std::cout << "done!" << std::endl;
	return 0;
}
