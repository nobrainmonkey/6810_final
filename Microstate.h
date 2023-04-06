// file: Microstate.h
//
// Header file for Microstate class
//
// Programmer: Xihe Han 
//
// Revision History:
//	04/05/2023 origional version
//
//***********************************


// ensure header file is only included once
#ifndef MICROSTATE_H
#define MICROSTATE_H

// include files
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include "Hamiltonian.h"
class Microstate{ 
	public:
		Microstate(int row, int col, double temp, hamiltonian_param_struct* hamiltonian_param);    //contructor

		~Microstate();    //destructor
						  
		// public variables
		double tempreature;
		int rows;
		int cols;
		hamiltonian_param_struct* hamiltonian_param_ptr;

		// access microstate matrix
		Eigen::MatrixXd* get_microstate_matrix_ptr();

		// print the microstate matrix
		void print_microstate_matrix();

		// graph the microstate matrix
		void graph_microstate_matrix();

		// graph the evolution of the microstat matrix
		void graph_evolve(int time, int dt);
	
		// evolve the microstate using MCMC seleciton rule:
		void evolve_microstate(int iteration);

	private:
		Eigen::MatrixXd* microstate_matrix_ptr; //microstate matrix that represents the 2-D lattice 
		
		// initialize the microstate matrix of given dimension row x col
		void initialize_microstate_matrix(Eigen::MatrixXd* microstate_matrix_ptr);
};

#endif
