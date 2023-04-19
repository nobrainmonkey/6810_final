// file: Microstate.h
//
// Header file for Microstate class
//
// Programmer: Xihe Han
//
// Revision History:
//	04/05/2023 original version
//  04/10/2023 added gradual termalize method
//  04/17/2023 changed parameters needed for graph_evolve
//***********************************

// ensure header file is only included once
#ifndef MICROSTATE_H
#define MICROSTATE_H

// include files
#include "Hamiltonian.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <omp.h>

class Microstate {
public:
  Microstate(int row, int col, double temp,
             hamiltonian_param_struct *hamiltonian_param); // contructor

  ~Microstate(); // destructor

  // public variables
  double temperature; // temperature that our microstate is thermalized to
  int rows;           // number of rows in the microstate lattice
  int cols;           // number of cols in the microstate lattice
  hamiltonian_param_struct *hamiltonian_param_ptr;

  // access microstate matrix
  Eigen::MatrixXd *get_microstate_matrix_ptr();

  // print the microstate matrix
  void print_microstate_matrix();

  // graph the microstate matrix
  void graph_microstate_matrix();

  // graph the evolution of the microstat matrix
  void graph_evolve(double time, int evolve_iteration, int fps);

  // evolve the microstate using MCMC seleciton rule:
  void evolve_microstate(int iteration);

  // evolve the microstate gradually with a start temp and a target temp
  void evolve_microstate_gradual(int iteration);

private:
  Eigen::MatrixXd *microstate_matrix_ptr; // microstate matrix that represents
                                          // the 2-D lattice

  // initialize the microstate matrix of given dimension row x col
  void initialize_microstate_matrix(Eigen::MatrixXd *microstate_matrix_ptr);

  // evolve the microstate once
  void evolve_microstate_once();
};

#endif
