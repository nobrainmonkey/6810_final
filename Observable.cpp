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

// return the ising energy of the single spin
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
		int num_average = 1000;
		// get the total amount of energy over the num_average iterations of the microstate.
		double total_energy = 0;
		for(int i=0; i < num_average; i++){
			microstate_ptr->evolve_microstate(1);
			total_energy += (Hamiltonian::hamiltonian_periodic_ising_microstate(microstate_ptr->get_microstate_matrix_ptr(), hamiltonian_param_ptr));
		}
		double microstate_total_energy = total_energy / double(num_average);
		double microstate_spin_energy = microstate_total_energy / double((row * col));
		total_sample_energy += microstate_spin_energy;
		delete microstate_ptr;
	}
	return total_sample_energy/double(sample_size);
}

// get the heat capacity by using dE/dT
// where the derivative is approximated by central derivate.
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
		int mc_steps = 1000;
		for(int i = 0; i < mc_steps; i++){
			microstate_ptr->evolve_microstate(1);
			double energy = Hamiltonian::hamiltonian_periodic_ising_microstate(microstate_ptr->get_microstate_matrix_ptr(), hamiltonian_param_ptr);
			E += energy;
			Esq += energy * energy;
		}
		double E_avg = E / mc_steps;
		double Esq_avg = Esq / mc_steps;

		double specific_heat = (Esq_avg - E_avg * E_avg) / (T * double(row * col)* hamiltonian_param_ptr->J);
		total_sample_heat_capacity += specific_heat;
		delete microstate_ptr;
	}
	return total_sample_heat_capacity / (double)sample_size;

}

double Observable::get_ising_entropy(double energy, double specific_heat, double T)
{
	double heat_capacity = energy/T + std::log(specific_heat);
	return heat_capacity;
}

int main()
{
	omp_set_num_threads(16);
	hamiltonian_param_struct hamiltonian_param;
	hamiltonian_param.J = 2.;
	hamiltonian_param.h = 0.0;
	int row = 20;
	int col = 20;
	int iteration = 4000000;
	int sample = iteration / 3;

	std::ofstream my_out;
	my_out.open("output.dat");
	my_out << "tempreature" << " " << std::setw(28) 
		   << "E"  << std::setw(37) << " "
		   << "Cv" << " " << std::setw(36)
		   << "S" << std::endl;
	for (double temp = 1.5; temp < 24.; temp += 0.1)
	{
		double energy = Observable::get_ising_energy(row, col, temp, iteration, sample, &hamiltonian_param);
		double cv = Observable::get_ising_heat_capacity(row, col, temp, iteration, sample, &hamiltonian_param);
		double entropy = Observable::get_ising_entropy(energy, cv, temp); 
		my_out <<std::scientific << std::setprecision(12)
			  << temp/hamiltonian_param.J << std::setw(20) <<" " 
			  << energy <<  std::setw(20) << " "
			  << cv  << std::setw(20) << " "
			  << entropy <<std::setw(20) << " "
			  <<std::endl; 
	}

	std::cout << "done!" << std::endl;
	return 0;
}
