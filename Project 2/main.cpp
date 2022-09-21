#include <iostream>
#include <cmath>
#include "utils.hpp"
#include "algorithm.hpp"

int main()
{
  
  int N=6; // number of equations to solve
  double h = 0.01; //stepsize
  double pi = 3.14159265358979323846; 

  // make tridiagonal matrix
  arma::mat A = arma::mat(N, N).fill(0.);
  A.diag() = arma::vec(N).fill(2./(h*h));
  A.diag(-1) = arma::vec(N-1).fill(-1./(h*h));
  A.diag(1) = arma::vec(N-1).fill(-1./(h*h));

  // solve for v and lambda using armadillo
  arma::vec eigval;
  arma::mat eigvec;

  arma::eig_sym(eigval, eigvec, A);
  eigval.print("Armadillo eigenvalues: ");

  arma::vec lambda(N);
  arma::mat V = arma::mat(N,N);

  for(int i = 1; i < N+1; i++)
  {
    lambda(i-1) = 2./(h*h) - 2/(h*h) * cos(i*pi / (N+1));

    for(int j=1; j < N+1; j++)
    {
      V(j-1, i-1) = sin(j*i*pi/(N+1));
    }
  }
  lambda.print("Analytical eigenvalues: ");

  // normalise to compare
  eigvec = arma::normalise(eigvec);
  eigvec.print("Armadillo eigenvector: ");
  V = arma::normalise(V);
  V.print("Analytic eigenvectors:");
 
  return 0;
}