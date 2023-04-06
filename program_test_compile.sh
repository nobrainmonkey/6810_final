#!/bin/bash
g++ -c Hamiltonian.cpp -o Hamiltonian.o
g++ -fopenmp -c Microstate.cpp -o Microstate.o
g++ -fopenmp -o Microstate main.cpp Hamiltonian.o Microstate.o

