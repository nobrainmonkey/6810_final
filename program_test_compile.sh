#!/bin/bash
g++ -c Hamiltonian.cpp -o Hamiltonian.o
g++ -fopenmp -c Microstate.cpp -o Microstate.o
g++ -fopenmp -c Observable.cpp -o Observable.o
g++ -fopenmp -O3 -o Observable main.cpp Hamiltonian.o Microstate.o Observable.o

