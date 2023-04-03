#include <random>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <armadillo>

using namespace arma;
using namespace std;

const int dimension = 400;

void Randomize_unit_matrix(mat* A);

int main()
{ 
	mat A(dimension, dimension, fill::randu);
	Randomize_unit_matrix(&A);
	return 0;
}

void Randomize_unit_matrix(mat* A)
{
	for(int i=0; i<dimension; i++)
	{
		for(int j=0; j<dimension; j++)
		{
			if((*A) (i,j) <= 0.5)
			{
				(*A) (i,j) = -1;
			}
			else{
				(*A) (i,j) = 1;
			}
		}
	}
}
