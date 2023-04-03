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

Because the problem serves to solve a 2-D MCMC problem, `eigen`'s lack of high-dimensional array methods wouldn't be an issue. If I wanted to solve for higher dimensional arrays, I could use `Blitz++`.

#### Numerical Methods:

For numerical methods such as numerical integration and derivatives, I have decided to use `GSL`. I think `GSL` is already very optimized for speed and it has all the numerical methods I wanted. 

#### Graphs:

For graphing macroscopic observables, I decide to try out `matplotlib-cpp`. This decision is not final, as plotting is not very integral to my detailed abstraction of my project. I might write an wrapper for `gnuplot` instead.

For graphing the microstates, I actually want to achieve this via some simple terminal `ASCII` graphing by myself. This should be rather simple as I will just be representing a matrix with entries `-1` or `1`. Using black and white is very sufficient for this purpose.


## Project Design Details
