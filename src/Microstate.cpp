// file: Microstate.cpp
//
// Source file for Microstate class
//
// Programmer: Xihe Han
//
// Revision History:
//  04/05/2023 origional version
//  04/10/2023 added gradual thermalize method
//  04/13/2023 fixed some pre-factors
//  04/17/2023 changed the parameters of graph_evolve
//				to be more readable
//
//	TODO: Precalcualte the Boltzmann factors
//        to reduce computation time.
//***********************************

// include files
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <random>
#include <omp.h>
#include "Hamiltonian.h"
#include "Microstate.h"
#include <iostream>
#include <chrono>
#include <thread>

// declearing random number generator using mt19937
std::mt19937 rng(std::random_device{}());
std::mt19937 gen(rng);

// construtor
Microstate::Microstate(int row, int col, double temp, hamiltonian_param_struct *hamiltonian_param)
{
	tempreature = temp;
	rows = row;
	cols = col;
	microstate_matrix_ptr = new Eigen::MatrixXd(rows, cols);
	initialize_microstate_matrix(microstate_matrix_ptr);
	hamiltonian_param_ptr = hamiltonian_param;
}

// destructpr
Microstate::~Microstate()
{
	delete microstate_matrix_ptr;
}

// initialize microstate matrix
// the matrix is initialized with dimension rows x cols
// with each element randomly assigned as -1. or 1.
void Microstate::initialize_microstate_matrix(Eigen::MatrixXd *microstate_matrix_ptr)
{
	std::bernoulli_distribution dis(0.5); // discrete probablity of 0.5 for 1. This means that we also have probablity 0.5 for -1
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			(*microstate_matrix_ptr)(i, j) = dis(gen) ? double(1.) : double(-1.); // shorthand for if-else statement
		}
	}
}

// getter function for the microstate matrix
Eigen::MatrixXd *Microstate::get_microstate_matrix_ptr()
{
	return microstate_matrix_ptr;
}

// print the microstate matrix
void Microstate::print_microstate_matrix()
{
	std::cout << "Microstate Matrix:" << std::endl;
	for (int j = 0; j < cols; j++)
	{
		std::cout << "| ";
		for (int i = 0; i < rows; i++)
		{
			std::cout << (*microstate_matrix_ptr)(i, j) << " ";
		}
		std::cout << "|" << std::endl;
	}
}

// graph the microstate matrix
// when the value of the matrix element is -1, we graph a space
// and if the value of the matrix element is 1, we graph a white square
// where the white sqaure is a unicode character

void Microstate::graph_microstate_matrix()
{
	std::cout << "Microstate Matrix:" << std::endl;
	for (int j = 0; j < cols; j++)
	{
		for (int i = 0; i < rows; i++)
		{
			if ((*microstate_matrix_ptr)(i, j) == -1.)
			{
				std::cout << "  "; // print space for -1
			}
			else
			{
				std::cout << "██"; // print white square for 1
			}
		}
		std::cout << std::endl;
	}
}

// routine to perform the Matropolis algorithm on a microstate.
// this routines only update a current microstate
void Microstate::evolve_microstate(int iteration)
{
	double inverseT = 1. / tempreature;
	// set random distrubution for row and col
	thread_local std::uniform_int_distribution<int> dis_row(0, rows - 1);
	thread_local std::uniform_int_distribution<int> dis_col(0, cols - 1);

	for (int i = 0; i < iteration; i++)
	{
		// flip a random element in our microstate
		int rand_row = dis_row(gen);
		int rand_col = dis_col(gen);
		// we assume to flip the spin at the random location
		(*microstate_matrix_ptr)(rand_row, rand_col) *= -1.;
		// delta E is just 2 times the flipped energy.
		// doing it this way so we dont need to calculate delta E by going through the entire matrix
		double deltaE = 2. * Hamiltonian::hamiltonian_periodic_ising_element(rand_row, rand_col, microstate_matrix_ptr, hamiltonian_param_ptr);

		// only reject the change of Delta E > 0 with the probablity e^(-deltaE/T)
		if (deltaE > 0)
		{
			std::uniform_real_distribution<double> dis_double(0.0, 1.0);
			double random_double = dis_double(gen); // random double is betwen 0 and 1
			double boltzman_factor = exp(-(deltaE * inverseT));
			if (random_double > boltzman_factor)
			{
				(*microstate_matrix_ptr)(rand_row, rand_col) *= -1.;
			}
		}
	}
}

// routine to evolve the microstate gradually from 5* target temperature
// this routine should always be used when termalizing the microstate
void Microstate::evolve_microstate_gradual(int iteration)
{
	double start_temp = 5 * tempreature; 
	double target_temp = tempreature;
	double inverseT = 1. / start_temp;
	double temp_step = (start_temp - target_temp) / iteration; //calculate the step of temperature given iteration 
	double current_temp = start_temp;
	int post_cooling_iteration = rows * cols * 100;    //number of iterations to run after target temp is reached.

	// Set random distribution for row and col
	thread_local std::uniform_int_distribution<int> dis_row(0, rows - 1);
	thread_local std::uniform_int_distribution<int> dis_col(0, cols - 1);

	// Evolve the microstate while cooling gradually
	for (int i = 0; i < iteration; i++)
	{
		// Flip a random element in our microstate
		int rand_row = dis_row(gen);
		int rand_col = dis_col(gen);
		(*microstate_matrix_ptr)(rand_row, rand_col) *= -1.;
		double deltaE = 2. * Hamiltonian::hamiltonian_periodic_ising_element(rand_row, rand_col, microstate_matrix_ptr, hamiltonian_param_ptr);

		// Only reject the change of Delta E > 0 with the probability e^(-deltaE/T)
		if (deltaE > 0)
		{
			std::uniform_real_distribution<double> dis_double(0.0, 1.0);
			double random_double = dis_double(gen);
			double boltzmann_factor = exp(-(deltaE * inverseT));
			if (random_double > boltzmann_factor)
			{
				(*microstate_matrix_ptr)(rand_row, rand_col) *= -1.;
			}
		}

		// Decrease temperature
		current_temp -= temp_step;
		inverseT = 1. / current_temp;
	}

	// Evolve the microstate at the target temperature for post_cooling_iteration iterations
	inverseT = 1. / target_temp;
	for (int i = 0; i < post_cooling_iteration; i++)
	{
		// Flip a random element in our microstate
		int rand_row = dis_row(gen);
		int rand_col = dis_col(gen);
		(*microstate_matrix_ptr)(rand_row, rand_col) *= -1.;
		double deltaE = 2. * Hamiltonian::hamiltonian_periodic_ising_element(rand_row, rand_col, microstate_matrix_ptr, hamiltonian_param_ptr);

		// Only reject the change of Delta E > 0 with the probability e^(-deltaE/T)
		if (deltaE > 0)
		{
			std::uniform_real_distribution<double> dis_double(0.0, 1.0);
			double random_double = dis_double(gen);
			double boltzmann_factor = exp(-(deltaE * inverseT));
			if (random_double > boltzmann_factor)
			{
				(*microstate_matrix_ptr)(rand_row, rand_col) *= -1.;
			}
		}
	}
}

// graph the evolution process
// this implementation is rather janky as I needed to use sleep to halt the program
// time is the total time in seconds where the graph will show up
// with a frame rate of field fps.
void Microstate::graph_evolve(double time, int evolve_iteration, int fps)
{    
	double t = 0.0;
	int dt_in_ms = 1000. / fps ;
	while (t < time * 1000)
	{
		// Update the microstate
		evolve_microstate(evolve_iteration);

		// Graph the microstate matrix
		graph_microstate_matrix();
		// Wait for some time
		std::this_thread::sleep_for(std::chrono::milliseconds(dt_in_ms));
		t += dt_in_ms;
		system("clear");
	}
}
