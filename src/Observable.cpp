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
#include "Observable.h"
#include <omp.h>
#include <iomanip>

// integrating using simspons rule if the number of points in x is odd
// and when the number of points is even, this routine is back to trapozodial rule
double integrate(std::vector<double> x, std::vector<double> y)
{
	int n = x.size();
	double h = (x[n - 1] - x[0]) / (n - 1);
	double sum_odd = 0.0, sum_even = 0.0;
	for (int i = 1; i < n - 1; i += 2)
	{
		sum_odd += y[i];
	}
	for (int i = 2; i < n - 1; i += 2)
	{
		sum_even += y[i];
	}
	double integral = (h / 3.0) * (y[0] + y[n - 1] + 4.0 * sum_odd + 2.0 * sum_even);
	return integral;
}
// return the ising energy of the single spin by doing the folowing:
// 1. initialize a microstate, evolve the microstate for `iteration` times so it reaches thermal equalibrium
// 2. keep evolving the microstate for `mc_steps` amount of time, record the energy for each evolution, and average the energy
// 3. return to step one for `sample_size` amount of time, record the averaged energy each time, and finally average all the samples
// 4. return the total sample energy per spin in the unit of J.
double Observable::get_ising_energy(int row, int col, double T, int iteration, int sample_size, hamiltonian_param_struct *hamiltonian_param_ptr)
{
	double total_sample_energy = 0.;
#pragma omp parallel for
	for (int i = 0; i < sample_size; i++)
	{
		Microstate *microstate_ptr = new Microstate(row, col, T, hamiltonian_param_ptr);

		// evolve the microstate to reach equalibruim
		microstate_ptr->evolve_microstate_gradual(iteration);

		// create a list to store energy after some iterations
		// we choose a list because it is dynamically sizedomp_set_num_threads(8);
		// and thread safe in case we ant to use omp.

		// after reaching equalibrium, we want to https://discordapp.com/channels/@me/1028655396090544259 further
		// evolve it 1000 time and record energy for each time
		// and take the final average.
		int mc_steps = iteration / 10; // amount of step to take average
									  // get the total amount of energy over the num_average iterations of the microstate.
		double total_energy = 0;
		for (int i = 0; i < mc_steps; i++)
		{
			microstate_ptr->evolve_microstate(1);
			total_energy += (Hamiltonian::hamiltonian_periodic_ising_microstate(microstate_ptr->get_microstate_matrix_ptr(), hamiltonian_param_ptr));
		}
		double microstate_total_energy = total_energy / double(mc_steps);
		double microstate_spin_energy = microstate_total_energy / double((row * col));
		total_sample_energy += microstate_spin_energy;
		delete microstate_ptr;
	}
	return total_sample_energy / (double(sample_size));
}

// return the ising heat capacity of the single spin by doing the folowing:
// 1. initialize a microstate, evolve the microstate for `iteration` times so it reaches thermal equalibrium
// 2. keep evolving the microstate for `mc_steps` amount of time, record the energy for each evolution, and use Cv = (<E^2> - <E>^2)/T to calcualte C_v
// 3. return to step one for `sample_size` amount of time, record the averaged energy each time, calculate Cv, and finally average all the samples
// 4. return the total sample specific heat per spin in the unit of J.
double Observable::get_ising_heat_capacity(int row, int col, double T, int iteration, int sample_size, hamiltonian_param_struct *hamiltonian_param_ptr)
{
	double total_sample_heat_capacity = 0;
#pragma omp parallel for
	for (int i = 0; i < sample_size; i++)
	{
		Microstate *microstate_ptr = new Microstate(row, col, T, hamiltonian_param_ptr);
		microstate_ptr->evolve_microstate_gradual(iteration); // bring system into thermal equalibrium
		double E = 0.;
		double Esq = 0.;
		int mc_steps = iteration / 10; // steps after thermal equalibrium
		for (int i = 0; i < mc_steps; i++)
		{
			microstate_ptr->evolve_microstate(1);
			double energy = Hamiltonian::hamiltonian_periodic_ising_microstate(microstate_ptr->get_microstate_matrix_ptr(), hamiltonian_param_ptr);
			E += energy;
			Esq += energy * energy;
		}
		double E_avg = E / mc_steps;
		double Esq_avg = Esq / mc_steps;
		double specific_heat = (Esq_avg - E_avg * E_avg) / (T * T * double(row * col));
		total_sample_heat_capacity += specific_heat;
		delete microstate_ptr;
	}
	return total_sample_heat_capacity / (double)sample_size;
}

// return the ising entropy per particle per J of the single spin by using the equation
// S = U/T + ln(Cv). Note the user should use the respective method to calculate E and Cv.
double Observable::get_ising_entropy(std::vector<double> T, std::vector<double> Cv)
{
	int size = T.size();
	std::vector<double> integrand; // integrand is Cv/T
								   // prepare the integrand with given Cv and T
	for (int i = 0; i < size; i++)
	{
		integrand.push_back(Cv[i] / T[i]);
	}
	double entropy = integrate(T, integrand);
	return entropy;
}

double Observable::get_ising_m(int row, int col, double T, int iteration, int sample_size, hamiltonian_param_struct *hamiltonian_param_ptr)
{
	double total_sample_m = 0.;
#pragma omp parallel for
	for (int i = 0; i < sample_size; i++)
	{
		Microstate *microstate_ptr = new Microstate(row, col, T, hamiltonian_param_ptr);

		// evolve the microstate to reach equalibruim
		microstate_ptr->evolve_microstate_gradual(iteration);

		// create a list to store energy after some iterations
		// we choose a list because it is dynamically sized
		// and thread safe in case we ant to use omp.

		// after reaching equalibrium, we want to https://discordapp.com/channels/@me/1028655396090544259 further
		// evolve it 1000 time and record energy for each time
		// and take the final average.
		int mc_steps = iteration / 10; // amount of step to take average
									  // get the total amount of energy over the num_average iterations of the microstate.
		double total_m = 0;
		for (int i = 0; i < mc_steps; i++)
		{
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
		Microstate *microstate_ptr = new Microstate(row, col, T, hamiltonian_param_ptr);
		microstate_ptr->evolve_microstate_gradual(iteration); // bring system into thermal equalibrium
		double m = 0.;
		double msq = 0.;
		int mc_steps = iteration / 10; // steps after thermal equalibrium
		for (int i = 0; i < mc_steps; i++)
		{
			microstate_ptr->evolve_microstate(1);
			double magnetization = microstate_ptr->get_microstate_matrix_ptr()->sum();
			m += magnetization;
			msq += magnetization * magnetization;
		}
		double m_avg = m / mc_steps;
		double msq_avg = msq / mc_steps;
		double chi = (msq_avg - m_avg * m_avg) / (double(row * col)*T);
		total_sample_chi += chi;
		delete microstate_ptr;
	}
	return total_sample_chi / (double)sample_size;
}
