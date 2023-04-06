// file: Microstate.cpp
//
// Source file for Microstate class
//
// Programmer: Xihe Han 
//
// Revision History:
//  04/05/2023 origional version
//
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
// construtor 
Microstate::Microstate(int row, int col, double temp, hamiltonian_param_struct* hamiltonian_param){
	tempreature = temp;
	rows = row;
	cols = col;
	microstate_matrix_ptr = new Eigen::MatrixXd(rows, cols);
	initialize_microstate_matrix(microstate_matrix_ptr);
	hamiltonian_param_ptr = hamiltonian_param;
}


//destructpr
Microstate::~Microstate() {
	delete microstate_matrix_ptr;
}

// initialize microstate matrix
// the matrix is initialized with dimension rows x cols
// with each element randomly assigned as -1. or 1.
void Microstate::initialize_microstate_matrix(Eigen::MatrixXd* microstate_matrix_ptr){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::bernoulli_distribution dis(0.5); //discrete probablity of 0.5 for 1. This means that we also have probablity 0.5 for -1
	#pragma omp parallel for
	for (int i=0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			(*microstate_matrix_ptr)(i,j) = dis(gen) ? 1 : -1; //shorthand for if-else statement 
		}
	}
}

// getter function for the microstate matrix
Eigen::MatrixXd* Microstate::get_microstate_matrix_ptr(){
	return microstate_matrix_ptr;
}

// print the microstate matrix
void Microstate::print_microstate_matrix(){
	std::cout << "Microstate Matrix:" << std::endl;
	for (int i = 0; i < rows; i++) {
		std::cout << "| ";
		for (int j = 0; j < cols; j++) {
			std::cout << (*microstate_matrix_ptr)(i, j) << " ";
		}
		std::cout << "|" << std::endl;
	}
}

// graph the microstate matrix
// when the value of the matrix element is -1, we graph a space
// and if the value of the matrix element is 1, we graph a white square
// where the white sqaure is a unicode character

void Microstate::graph_microstate_matrix(){
	std::cout << "Microstate Matrix:" << std::endl;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if ((*microstate_matrix_ptr)(i, j) == -1) {
				std::cout << "  "; // print space for -1
			} else {
				std::cout << "██"; // print white square for 1
			}
		}
		std::cout << std::endl;
	}
}

// code to evolve the microstate with given iteration using single step monte carlo markov chain
// openmp should be implemented here
void Microstate::evolve_microstate(int iteration){
	double inverseT = 1./tempreature; 
	//set random distrubution for row and col
	thread_local std::mt19937 gen(std::random_device{}());
	thread_local std::uniform_int_distribution<int> dis_row(0, rows-1);
	thread_local std::uniform_int_distribution<int> dis_col(0, cols-1);

	#pragma omp parallel for
	for (int i =0; i < iteration; i ++)
	{
		//flip a random element in our microstate
		int rand_row = dis_row(gen);
		int rand_col = dis_col(gen);
		(*microstate_matrix_ptr)(rand_row, rand_col) *= -1;
		double deltaE = 2. * Hamiltonian::hamiltonian_periodic_ising_element(rand_row, rand_col, microstate_matrix_ptr, hamiltonian_param_ptr);

		//only reject the change of Delta E > 0 with the probablity e^(-deltaE/T)
		if(deltaE > 0)
		{
			std::uniform_real_distribution<double> dis_double(0.0, 1.0);
			double random_double = dis_double(gen);		//random double is betwen 0 and 1
			double boltzman_factor = exp(-(deltaE* inverseT));
			if(random_double > boltzman_factor)
			{
				(*microstate_matrix_ptr)(rand_row, rand_col) *= -1;
			}
		}

	}
}

// graph the evolution process
// this implementation is rather janky as I needed to use sleep to halt the program
// The total amount of iteration this will graph is time / dt -1
void Microstate::graph_evolve(int time, int dt) {
    double t = 0.0;
    while (t < time) {
        // Update the microstate
        evolve_microstate(50);

        // Graph the microstate matrix
        graph_microstate_matrix();
        // Wait for some time
		std::this_thread::sleep_for(std::chrono::milliseconds(dt));
        t += dt;
    }
}


int main()
{
	hamiltonian_param_struct hamiltonian_variables;
	hamiltonian_variables.J = 2;
	hamiltonian_variables.h = 0.;
	omp_set_num_threads(16);
	Microstate microstate = Microstate(20, 20, 10, &hamiltonian_variables);
	microstate.graph_microstate_matrix();
	microstate.evolve_microstate(1000000);
	microstate.graph_microstate_matrix();
	std::cout<<"done!";
}
