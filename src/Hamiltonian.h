//  file Hamiltonian.h
//
//  Header file for the Hamiltonian functions
//
//  Programmer: Xihe Han
//
//  Revision History:
//		04/05/2023 origional version
//
//***********************************************

// ensuring the header file is only defined once
#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

// defining the structure of the hamiltonian parameters.
// In the case of Ising model, all we need is J and h.
struct hamiltonian_param_struct
{
	double J;
	double h;
};

// include files
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>

// namespace definition and function prototypes
namespace Hamiltonian
{
	// Define the Hamiltonian function for a single spin in 2-D ising model
	// This function return the interaction energy of an element at (i-1, j-1)
	// with its neighbors with periodic boundary conditon.
	double hamiltonian_periodic_ising_element(int i, int j, Eigen::MatrixXd *microstate_matrix_ptr, hamiltonian_param_struct *hamiltonian_param);

	// Define the Hamiltonian fucntion for the entire microstate of a system of sice
	// NxM.
	double hamiltonian_periodic_ising_microstate(Eigen::MatrixXd *microstate_matrix_ptr, hamiltonian_param_struct *hamiltonian_param);

};

#endif
