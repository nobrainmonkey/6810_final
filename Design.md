# P6810 Final- Ising Model Simulation

Project GitHub Page: https://github.com/nobrainmonkey/6810_final

Author: Xihe Han

---

## Project Objectives

### Tools Related Learning

- [X] Learn How to use a version control software connected with a remote repository.

### Computationally solve the problem of a 2-D two-level Ising model with a given Hamiltonian under a magnetic field $h$.

- [X] Microstate generation.
- [ ] Micro-sate visualization. (Maybe omitted)
- [ ] Micro-sate enumeration using Metropolis–Hastings algorithm (MCMC).
- [ ] Observables calculation.
- [ ] Parallel computing using `Openmp`.

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

The ultimate goal of this project is to **solve a physics problem in the best way I could think of.** This means that instead of re-inventing the wheel and writing my own `matrix` class and my own `numerical` packages, I do want to mainly use packages avaliable online for the best performance/accuracy purposes. I will probably implement some hand-written integration/derivative methods, but those are not going to be included in `main.cpp`- they are just to show that I learned some computational method in this course.

### Technology stack

#### Matrix Computation:

For matrix related computation, I have decided to use `eigen`. After doing some research on the internet, I learned that `eigen` is the most canonical choice for doing dense matrix operation and it is in many cases blazingly fast.

Because the problem serves to solve a 2-D MCMC problem, `eigen`'s lack of high-dimensional array methods wouldn't be an issue. If I wanted to solve for higher dimensional arrays, I could use `blitz++`.

#### Numerical Methods:

For numerical methods such as numerical integration and derivatives, I have decided to use `GSL`. I think `GSL` is already very optimized for speed and it has all the numerical methods I wanted.

#### Graphs:

For graphing macroscopic observables, I decide to try out `matplotlib-cpp`. This decision is not final, as plotting is not very integral to my detailed abstraction of my project. I might write an wrapper for `gnuplot` instead. My PI will make me learn `root` eventually, so that is probably also an option.

For graphing the microstates, I actually want to achieve this via some simple terminal `ASCII` graphing by myself. This should be rather simple as I will just be representing a matrix with entries `-1` or `1`. Using black and white is very sufficient for this purpose.

## Project Design Details

### Object Oriented Programming

Because we are using `c++`, we should embrace the modularity of `OOP`. This means that we seek to maximize the readability and adaptability of our code. Below are the necessary objects and their structure for this program.

#### Microstate

The `Microstate` class should inlucde everything we need to describe the microstate of out spin system. This should include **a complete microscopic description of our microstate.** However, any macroscopic quantities will Not be covered here. Because we are dealing with a 2-D lattice, we shall use a 2-D array from `eigen` to represent this state.

```

```


#### Hamiltonian

The main function of our Hamiltonian class is to **calcuate the energy of a given microsate.** Because energy is a macroscopic quantity, we probably will never deriectly utilize this class in our `main` function, but this Hamiltonian will be used for our macroscopic quantities' calculation.

```
class Hamiltonian (microstate* micro_state)
private: double J, doubhe h
double J = micro_state.get(J)
doubel h = micro_state.get(h)


double Get_element_energy(microstate* microstate, int i, int j)
// return the ith row, jth col spin energy calculation using E = - J sum_<ij> S_i S_j - h* S_i
// where <ij> stands for the nearest neightbor of spin located at (i,j) with periodic boundary condition

double Get_microstate_energy(microstate* microstate)
//return the total energy of the entire microstate by intrearting through every element of the microstate 
//and calling Get_element_energy()

```

Note: 

* The nearest neightbor of periodic B.C. will probably be some what difficult. Now all I can think of is using `mod` operator. I might seperate this part out to another function in this class.



* We start from a random microstate.
* We flip a random spin in that microstate.
* If the microstate energy becomes lower, we accept this change.
* If the microstate energy becomes higher, we accept this change with probablity $e^{-\Delta E/T }$
* Iterate for a large number until we reach equalibrium.

* We start from a random microstate.
* We flip a random spin in that microstate.
* If the microstate energy becomes lower, we accept this change.
* If the microstate energy becomes higher, we accept this change with probablity $e^{-\Delta E/T }$
* Iterate for a large number until we reach equalibrium.
