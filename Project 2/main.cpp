#include <iostream>
#include <cmath>
#include "utils.hpp"
#include "algorithm.hpp"

int main()
{

  int N=6; // number of equations to solve
  double h = 0.01; //stepsize
  double d = 2./(h*h); // main diagonal
  double a = -1./(h*h); // upper + lower diagonal

  // make tridiagonal matrix
  arma::mat A = make_tridiag(N, a, d);

  // solve for v and lambda using armadillo
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, A);
  eigvec = arma::normalise(eigvec); // normalise for comparison with analytic solution
  
  // find analytic solution
  arma::vec lambda(N); //eigenvalues
  arma::mat V = arma::mat(N,N); //eigenvectors
  solve_analytic(N, a, d, lambda, V);

  // now we compare the solutions
  arma::vec are_equal(N); 
  bool eigvals_equal = arma::approx_equal(lambda, eigval, "absdiff", 1e-8);

  for(int i=0; i<N; i++)
  {
    are_equal(i) = arma::approx_equal(eigvec.col(i), V.col(i), "absdiff", 1e-8) || arma::approx_equal(eigvec.col(i), -1*V.col(i), "absdiff", 1e-8);
  }
  bool eigvecs_equal = arma::all(are_equal);

  if(eigvals_equal && eigvecs_equal)
  {
    std::cout << "The eigenvalues and eigenvalues are equal." << std::endl;
  } else{
    std::cout << "The eigenvalues and aigenvalues are NOT equal." << std::endl;
  }
 
  return 0;
}