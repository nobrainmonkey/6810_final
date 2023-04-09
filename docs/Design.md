# P6810 Final- Ising Model Simulation

Project GitHub Page: https://github.com/nobrainmonkey/6810_final

Author: Xihe Han

---

## Project Objectives

### Tools Related Learning

- [X] Learn How to use a version control software connected with a remote repository.

### Computationally solve the problem of a 2-D two-level Ising model with a given Hamiltonian under a magnetic field $h$.

- [X] Hamiltonian construction.
- [X] Microstate generation.
- [X] Micro-sate visualization. (Maybe omitted)
- [X] Micro-sate enumeration using Metropolis–Hastings algorithm (MCMC).
- [ ] Observables calculation.
- [X] Parallel computing using `Openmp`.

### Graphical Representation

- [ ] Graph calculated observables under different conditions.
- [ ] Maybe explore more graphing software than `gnuplot`. As my research area is in nuclear physics, maybe I will try out `root`.

### Error Analysis

- [ ] Microstate operations error analysis.
- [ ] Observable calculation error analysis.
- [ ] Total error analysis.

### Writing

- [ ] Summary, write-up, and potential improvements.

---

## Project Design Overview

### Design Philosophy

The ultimate goal of this project is to **solve a physics problem in the best way I could think of.** This means that instead of re-inventing the wheel and writing my own `matrix` class and my own `numerical` packages, I do want to mainly use packages available online for the best performance/accuracy purposes. I will probably implement some handwritten integration/derivative methods, but those are not going to be included in `main.cpp`- they are just to show that I learned some computational method in this course.

### Technology stack

#### Matrix Computation:

For matrix related computation, I have decided to use `eigen`. After doing some research on the internet, I learned that `eigen` is the most canonical choice for doing dense matrix operation, and it is in many cases fast.

Because the problem serves to solve a 2-D MCMC problem, `eigen`'s lack of high-dimensional array methods wouldn't be an issue. If I wanted to solve for higher dimensional arrays, I could use `blitz++`.

#### Numerical Methods:

For numerical methods such as numerical integration and derivatives, I have decided to use `GSL`. I think `GSL` is already very optimized for speed ,and it has all the numerical methods I wanted.

#### Graphs:

For graphing macroscopic observables, I decide to try out `matplotlib-cpp`. This decision is not final, as plotting is not very integral to my detailed abstraction of my project. I might write a wrapper for `gnuplot` instead. My PI will make me learn `root` eventually, so that is probably also an option.

For graphing the microstates, I actually want to achieve this via some simple terminal `ASCII` graphing by myself. This should be rather simple as I will just be representing a matrix with entries `-1` or `1`. Using black and white is very sufficient for this purpose.

## Project Design Details

### Object-Oriented Programming

Because we are using `c++`, we should embrace the modularity of `OOP`. This means that we seek to maximize the readability and adaptability of our code. Below are the necessary objects and their structure for this program.

#### Microstate

The `Microstate` class should include everything we need to describe the microstate of our spin system. This should include **a complete microscopic description of our microstate.** However, any macroscopic quantities will Not be covered here. Because we are dealing with a 2-D lattice, we shall use a 2-D array from `eigen` to represent this state.  To define this lattice, we must define the row number and col number.

The `Microstate` should also be able to evolve according to the Metropolis algorithm:

1. We start from a random microstate.
2. We flip a random spin in that microstate.
3. If the microstate energy becomes lower, we accept this change.
4. If the microstate energy becomes higher, we accept this change with probability $e^{-\Delta E/T }$
5. Iterate for a large number until we reach equilibrium.

The `Microstate` class should have the following structure.

```cpp
#include 'Haimtonian.h'
#include 'iostream'
#include 'cmath'
#include 'random'
#include <omp.h>

class Microstate (int row, int col, double temp)


// Hamiltonian class need to access the number of row and col.
private: 
int rownum = row 
int colnum = col 
double T = temp
double inverseT = 1. / temp
omp_set_num_threads(16);


// function prototype for initialize microsate
mat* Initialize_matrix(int row, int col)

// Create initial matrix 
mat* microstate_matrix = new mat(rownum, colnum, randomly fill between -1 and 1)
Initialize_matrix(microstate_matrix)

// getter and setter function for row, col, and microstate_matrix




// Initialize the microsotate to be a row X col matrix with random values from -1 to 1 defined by eigen
void Initialize_matrix(mat* A){
for (i=0; i < row; i++){
	for(j = 0; j < col; j++){
		if A(i, j) < = 0 -> A(i,j) = -1
		else A(i,j) = 1
	}
}
return * mat A
}



// code to evolve the microstate with given iteration using single step monte carlo markov chain
// openmp should be implemented here

void evolve_microstate(int iteration){
 
	//set random distrubution for row and col
	thread_local std::mt19937 gen(std::random_device{}());
    	thread_local std::uniform_int_distribution<int> dis_row(0, row_num-1);
   	thread_local std::uniform_int_distribution<int> dis_col(0, col_num-1);

	#pragma omp parallel for
	for (int i =0; i < iteration; i ++)
	{
		//flip a random element in our microstate
		int rand_row = dis_row(gen)
		int rand_col = dis_col(gen)
		*microstate_matrix(rand_row, rand_col) *= -1;
		deltaE = 2. * Get_ising_element_energy(rand_row, rand_col, *microstate_matrix)

		//only reject the change of Delta E > 0 with the probablity e^(-deltaE/T)
		if(deltaE > 0)
		{
			std::uniform_real_distribution<double> dis_double(0.0, 1.0);
			double random_double = dis_double(gen);		//random double is betwen 0 and 1
			double boltzman_factor = exp(-(deltaE* inverseT)
			if(random_double > boltzman_factor)
			{
				*microstate_matrix(rand_row, rand_col) *= -1;
			}
		}

	}
}

// graphing function should go here 
void print_microstate_matrix()
{
	// go through the matrix, print each element in a nice way 
}


```

#### Hamiltonian

The main function of our Hamiltonian class is to **calculate the energy of a given microsate.** Because energy is a macroscopic quantity, we probably will never directly utilize this class in our `main` function, but this Hamiltonian will be used for our macroscopic quantities' calculation. This probably won't be a class but merely a header file and source file for other classes.

```cpp
Hamiltonian.h

struct hamiltonian_param_struct{
	// Parameters regarding H should go here. For Ising model case, we have 
	double J;
	double h;
}

double Get_ising_element_energy(mat* microstate_matrix, int i, int j, hamiltonian_param_struct* hamiltonian_param)
// return the ith row, jth col spin energy calculation using hamiltonian_ising_function
// where <ij> stands for the nearest neightbor of spin located at (i,j) with periodic boundary condition

double Get_ising_microstate_energy(mat* microstate_matrix, hamiltonian_param_struct* hamiltonian_param)
//return the total energy of the entire microstate by intrearting through every element of the microstate 
//and calling Get_element_energy()


```

Note:

* The nearest neighbors of periodic B.C. will probably be some-what difficult. Now all I can think of is using `mod` operator. I might separate this part out to another function in this class.

#### Macro_Quantity

This is nothing but a collection of function that calculates the macroscopic quantities. It should have the following structure:

```cpp
#include <Hamiltonian.h>
#include <gsl/>

// get energy of the microstate 
double get_Energy(Microstate * microstate) {
	return Hamiltonian::Get_ising_microstate_energy(microstate->microstate_matrix);
}

// get specific head of mirostate over a range of temp 
double get_specific_heat(double T, double deltaT, int iteration){
C = Hamiltonian:: Get_ising_microstate_energy((new Microstate(N,M,T+deltaT).evole(iteration))->microstate_matrix) - Hamiltonian:: Get_ising_microstate_energy((new Microstate(N,M,T).evole(iteration))->microstate_matrix)/(deltaT)
return C
}

//calculate entropy
double get_entropy(double upper_T, int iteration){
return numerical_integral(CV(T)/T,T,0, upper_T)
}

//calcuate magnetization 
double get_magnetization(double T)
return total_Si

// calculate chi
double get_chi(double T, double deltaT)
return get_magnetization(T+ deltaT) - get_magnetization(T) / delta T


```

#### Main

Our main program should do the following:

1. Create a microstate.
2. Evolve the microstate.
3. (Optional) Graph the microstate.
4. Define a temperature mesh `T_mesh.`
5. Calculate the macroscopic quantities over `T_mesh.`
6. Graph the macroscopic quantities.

Thing to define in my main function

```cpp
hamiltonian_param_strcut // we define J and h as a field here
omp_set_num_threads // set number of threads to use for omp
row number, col number, Tempreature //for us to generate microstate of a given size and tempreature


// Maybe make a command line interface for us to change various settings!
```

# End
