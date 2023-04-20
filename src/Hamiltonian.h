//  file Hamiltonian.h
//
//  Header file for the Hamiltonian functions
//
//  Programmer: Xihe Han
//
//  Revision History:
//              04/05/2023 original version
//
//***********************************************

// ensuring the header file is only defined once
#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

// defining the structure of the Hamiltonian parameters.
// In the case of Ising model, all we need is J and h.
struct hamiltonian_param_struct {
  double J; // spin-spin interaction energy
  double h; // spin-field interaction energy
};

// include files
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>

// namespace definition and function prototypes
namespace Hamiltonian {
// Define the Hamiltonian function for a single spin in 2-D Ising model
// This function return the interaction energy of an element at (i-1, j-1)
// with its neighbors with periodic boundary condition.
double
hamiltonian_periodic_ising_element(int i, int j,
                                   Eigen::MatrixXd *microstate_matrix_ptr,
                                   hamiltonian_param_struct *hamiltonian_param);

// Define the Hamiltonian function for the entire microstate of a system of size
// NxM.
double hamiltonian_periodic_ising_microstate(
    Eigen::MatrixXd *microstate_matrix_ptr,
    hamiltonian_param_struct *hamiltonian_param);

}; // namespace Hamiltonian

#endif
