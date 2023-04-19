// file: Hamiltonian.cpp
//
// Definition for the Hamiltonian namespace functions
//
// Programmer: Xihe Han
//
// Revision history:
//	04/05/2023 Origianal version
//*****************************************************

// include files
#include "Hamiltonian.h"

// function definition in Hamiltonian.h
// return the energy calculation of the nearest neighbor interaction of a spin
// located at (i,j)
double Hamiltonian::hamiltonian_periodic_ising_element(
    int i, int j, Eigen::MatrixXd *microstate_matrix_ptr,
    hamiltonian_param_struct *hamiltonian_param) {
  double J = hamiltonian_param->J; // spin-spin interaction constant
  double h = hamiltonian_param->h; // magnetic field value

  int const rows = microstate_matrix_ptr->rows();
  int const cols = microstate_matrix_ptr->cols();

  // calculate the change in row index for up/down neighbor
  // and the change in col index for left/right neighbor
  int up = (i - 1 + rows) % rows;
  int down = (i + 1) % rows;
  int left = (j - 1 + cols) % cols;
  int right = (j + 1) % cols;

  // apply the equation E_ij = sum<S_ij> - h*S_i.
  // this equation comes from the hamiltonian of ising model
  double element_interaction_energy =
      -J * ((*microstate_matrix_ptr)(i, j) * (*microstate_matrix_ptr)(up, j) +
            (*microstate_matrix_ptr)(i, j) * (*microstate_matrix_ptr)(down, j) +
            (*microstate_matrix_ptr)(i, j) * (*microstate_matrix_ptr)(i, left) +
            (*microstate_matrix_ptr)(i, j) *
                (*microstate_matrix_ptr)(i, right)) -
      h * (*microstate_matrix_ptr)(i, j);

  return element_interaction_energy;
}

// function definition in Hamiltonian.h
// return the energy calculation for the entire microstate.
double Hamiltonian::hamiltonian_periodic_ising_microstate(
    Eigen::MatrixXd *microstate_matrix_ptr,
    hamiltonian_param_struct *hamiltonian_param) {
  int const rows = microstate_matrix_ptr->rows();
  int const cols = microstate_matrix_ptr->cols();

  double microstate_energy = 0.;

  // go through the microstate
  // and sum the energy for each spin.
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      microstate_energy += hamiltonian_periodic_ising_element(
          i, j, microstate_matrix_ptr, hamiltonian_param);
    }
  }

  return microstate_energy / 2.; // divide by 2 to avoid double counting.
}
