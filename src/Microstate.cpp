// file: Microstate.cpp
//
// Source file for Microstate class
//
// Programmer: Xihe Han
//
// Revision History:
//  04/05/2023 original version
//  04/10/2023 added gradual thermalize method
//  04/13/2023 fixed some pre-factors
//  04/17/2023 changed the parameters of graph_evolve
//                              to be more readable
//      04/18/2023 fixed multiple rng, made MC step modular
//
//      TODO: Pre-calculate the Boltzmann factors
//        to reduce computation time.
//***********************************

// include files
#include "Microstate.h"
#include "Hamiltonian.h"
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <iostream>
#include <omp.h>
#include <random>
#include <thread>

//***************************************************************************************************************************************
//**********************************************************RNG
// routines*****************************************************************
//***************************************************************************************************************************************

// declearing random number generator using mt19937
thread_local std::mt19937 global_rng(std::random_device{}

                                     ());

// generate a random integer between 0 and N-1
int random_integer(int N) {
  thread_local std::uniform_int_distribution<int> dist(0, N - 1);
  return dist(global_rng);
}

// generate a random double between 0 and 1
double random_unitary_double() {
  thread_local std::uniform_real_distribution<double> dist(0, 1);
  return dist(global_rng);
}

int random_sign() {
  thread_local std::uniform_int_distribution<int> dist(0, 1);
  int sign = dist(global_rng) * 2 - 1;
  return double(sign);
}

//***************************************************************************************************************************************
//******************************************************Class
// Initialization*************************************************************
//***************************************************************************************************************************************

// construtor
Microstate::Microstate(int row, int col, double temp,
                       hamiltonian_param_struct *hamiltonian_param) {
  temperature = temp;
  rows = row;
  cols = col;
  microstate_matrix_ptr = new Eigen::MatrixXd(rows, cols);
  initialize_microstate_matrix(microstate_matrix_ptr);
  hamiltonian_param_ptr = hamiltonian_param;
}

// destructpr
Microstate::~Microstate() { delete microstate_matrix_ptr; }

// initialize microstate matrix
// the matrix is initialized with dimension rows x cols
// with each element randomly assigned as -1. or 1.
void Microstate::initialize_microstate_matrix(
    Eigen::MatrixXd *microstate_matrix_ptr) {
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      (*microstate_matrix_ptr)(i, j) =
          double(random_sign()); // fill the matrix with random 1 or -1
    }
  }
}

// getter function for the microstate matrix
Eigen::MatrixXd *Microstate::get_microstate_matrix_ptr() {
  return microstate_matrix_ptr;
}

//***************************************************************************************************************************************
//******************************************************Graphing and
// Printing************************************************************
//***************************************************************************************************************************************

// print the microstate matrix
void Microstate::print_microstate_matrix() {
  std::cout << "Microstate Matrix:" << std::endl;
  for (int j = 0; j < cols; j++) {
    std::cout << "| ";
    for (int i = 0; i < rows; i++) {
      std::cout << (*microstate_matrix_ptr)(i, j) << " ";
    }
    std::cout << "|" << std::endl;
  }
}

// graph the microstate matrix
// when the value of the matrix element is -1, we graph a space
// and if the value of the matrix element is 1, we graph a white square
// where the white sqaure is a unicode character

void Microstate::graph_microstate_matrix() {
  std::cout << "Microstate Matrix:" << std::endl;
  for (int j = 0; j < cols; j++) {
    for (int i = 0; i < rows; i++) {
      if ((*microstate_matrix_ptr)(i, j) == -1.) {
        std::cout << "  "; // print space for -1
      } else {
        std::cout << "██"; // print white square for 1
      }
    }
    std::cout << std::endl;
  }
}

// graph the evolution process
// this implementation is rather janky as I needed to use sleep to halt the
// program time is the total time in seconds where the graph will show up with a
// frame rate of field fps.

void Microstate::graph_evolve(double time, int evolve_iteration, int fps) {
  double t = 0.0;
  int dt_in_ms = 1000. / fps;
  while (t < time * 1000) {
    // Update the microstate
    evolve_microstate(evolve_iteration);

    // Graph the microstate matrix
    graph_microstate_matrix();
    // Wait for some time
    std::this_thread::sleep_for(std::chrono::milliseconds(dt_in_ms));
    t += dt_in_ms;
    std::cout << "\033[2J\033[1;1H" << std::endl;
  }
}

//***************************************************************************************************************************************
//*********************************************************MCMC
// Algorithm****************************************************************
//***************************************************************************************************************************************

// source code to evolve microstate once
void Microstate::evolve_microstate_once() {
  int rand_row = random_integer(rows);
  int rand_col = random_integer(cols);

  (*microstate_matrix_ptr)(rand_row, rand_col) *= -1.;
  // delta E is just 2 times the flipped energy.
  // doing it this way so we dont need to calculate delta E by going through the
  // entire matrix
  double deltaE = 2. * Hamiltonian::hamiltonian_periodic_ising_element(
                           rand_row, rand_col, microstate_matrix_ptr,
                           hamiltonian_param_ptr);

  // only reject the change of Delta E > 0 with the probablity e^(-deltaE/T)
  if (deltaE > 0) {
    double boltzman_factor = std::exp(-(deltaE) / temperature);
    double random_double = random_unitary_double();
    if (random_double >
        boltzman_factor) // if the random double is greater than the boltzman
                         // factor, we flip spin back.
    {
      (*microstate_matrix_ptr)(rand_row, rand_col) *= -1.;
    }
  }
}

// routine to perform the Matropolis algorithm on a microstate for a given
// iteration this routines only update a current microstate
void Microstate::evolve_microstate(int iteration) {
#pragma omp parallel for
  for (int i = 0; i < iteration; i++) {
    evolve_microstate_once();
  }
}

// routine to evolve the microstate gradually from 5* target temperature
// this routine should always be used when termalizing the microstate
void Microstate::evolve_microstate_gradual(int iteration) {

  double target_temp = temperature;
  double current_temp = 5 * temperature;

  temperature = current_temp;

  double const TEMPERATURE_STEPS = 31;
  int const COOLING_ITERATION = rows * cols * 50;

  for (int i = 0; i < TEMPERATURE_STEPS; i++) {
    evolve_microstate(COOLING_ITERATION);
    current_temp *= 0.95;
    temperature = current_temp;
  }
  // Evolve the microstate at the target temperature for post_cooling_iteration
  // iterations
  temperature = target_temp;
  evolve_microstate(iteration);
}
