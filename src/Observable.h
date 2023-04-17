// Fila name: Observable.h
//
// Headerfile for calculating observables from a microstate
//
// Programmer: Xihe Han
//
// Version history:
//	04/06/2023: original version
//	04/13/2023: Added Cv, chi calculation
//
//**************************************************************

#ifndef OBSERVABLE_H
#define OBSERVABLE_H

// include files
#include "Microstate.h"
#include <vector>
// namespace definition and function prototypes
namespace Observable
{
	// function prototype for calculating
	// energy of a microstate of a given
	// dimension row x col and a given
	// hamiltonian parameters.
	double get_ising_energy(int row, int col, double T, int iteration, int sample_size, hamiltonian_param_struct *hamiltonian_param_ptr);

	// function prototype for calculating
	// heat capcity of a microstate of a given
	// dimension row x col by using numerical differentiation
	// dE/dT with step dT and a given hamiltonian parameters.
	double get_ising_heat_capacity(int row, int col, double T, int iteration, int sample_size, hamiltonian_param_struct *hamiltonian_param_ptr);

	// function prototype for calculating
	// entropy of a microstate of a given
	// dimension row x col by using numerical integration of Cv/T
	// with step dT and a given hamiltonian parameters
	double get_ising_entropy(std::vector<double> T, std::vector<double> Cv);

	// function prototype for calculating
	// magnetization of a microstate of a given
	// dimension row x col and a given hamiltonian parameters.
	double get_ising_m(int row, int col, double T, int iteration, int sample_size, hamiltonian_param_struct *hamiltonian_param_ptr);

	// function prototype for calculating
	// magnetization susceptibility of a microstate of a given
	// dimension row x col and a given hamiltonian parameters with step dh
	double get_ising_chi(int row, int col, double T, int iteration, int sample_size, hamiltonian_param_struct *hamiltonian_param_ptr);

}
#endif
