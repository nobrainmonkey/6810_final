#!/bin/bash
g++ -fopenmp -O3 -c Hamiltonian.cpp -o Hamiltonian.o
g++ -fopenmp -O3 -c Microstate.cpp -o Microstate.o
g++ -fopenmp -O3 -c Observable.cpp -o Observable.o
g++ -fopenmp -O3 -o Observable main.cpp Hamiltonian.o Microstate.o Observable.o

