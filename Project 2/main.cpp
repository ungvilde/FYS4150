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
  arma::vec are_equal(N); // vector containing boolean value when comparing eigenvectors
  bool eigvals_equal = arma::approx_equal(lambda, eigval, "absdiff", 1e-8);

  for(int i=0; i<N; i++)
  {
    // check if eigenvectors are the same (also checks vector when scaled by -1)
    are_equal(i) = arma::approx_equal(eigvec.col(i), V.col(i), "absdiff", 1e-8) || arma::approx_equal(eigvec.col(i), -1*V.col(i), "absdiff", 1e-8);
  }
  bool eigvecs_equal = arma::all(are_equal);

  if(eigvals_equal && eigvecs_equal)
  {
    std::cout << "The eigenvalues and eigenvalues are equal." << std::endl;
  } else{
    std::cout << "The eigenvalues and eigenvalues are NOT equal." << std::endl;
  }
 
  arma::mat A_test = arma::mat(4,4).fill(0.);
  A_test.diag() = arma::vec(4).fill(1.);
  A_test(3,0) = 0.5;
  A_test(2,1) = -0.7;
  A_test += A_test.t();
  A_test.print("A_test:"); 
  int l = 0;
  int k = 0;
  double max_elem = max_offdiag_symmetric(A_test, k, l);
  std::cout << "Max off-diagonal element: " << max_elem << std::endl;
  std::cout << "Indeces: " << k << ", " << l << std::endl;

  return 0;
}