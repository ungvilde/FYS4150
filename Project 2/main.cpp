#include <iostream>
#include <cmath>
#include "utils.hpp"
#include "algorithm.hpp"

int main()
{
  // Problem 2: Compare armadillo and analyitcal solution

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
  lambda.print("Analytic eigenvalues:");
  V.print("analytic eigenvectors:");

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

  // Problem 3: test max_offdiag_symmetric-function for plroblem 3
  arma::mat A_test = arma::mat(4,4).fill(0.);
  A_test.diag() = arma::vec(4).fill(1.);
  A_test(3,0) = 0.5;
  A_test(2,1) = -0.7;
  A_test += A_test.t();
  A_test.print("A_test:"); 
  int l = 0;
  int k = 0;
  double max_elem = max_offdiag_symmetric(A_test, k, l);
  std::cout << "Max. off-diagonal element: " << max_elem << std::endl;
  std::cout << "Indeces: " << k << ", " << l << std::endl;


  // Problem 4: test jacobi rotation algorithm
  // arma::vec eigenvalues;
  // arma::mat eigenvectors;
  int iterations;
  bool converged;

  A = make_tridiag(6, a, d);
  A.print("Original A:");
  //jacobi_eigensolver(A, 1e-4, eigenvalues, eigenvectors, 100, iterations, converged);
  int maxiter = 10000000;
  double eps = 1e-8;

  k = 0;
  l = 1;
  arma::mat R = arma::mat(N, N, arma::fill::eye);

  max_elem = max_offdiag_symmetric(A, k, l);

  //while(max_elem > eps)
  for(int i=0; i<=maxiter; i++)
  {
      if(max_elem < eps) // if convergence was reached before hitting maxiter
      {
          converged = true;
          iterations = i;
      }

      jacobi_rotate(A, R, k, l);
      //std::cout <<"k=" << k << ", l ="<< l << std::endl;
      max_elem = max_offdiag_symmetric(A, k, l);
      
      if(i==maxiter)
      {
          iterations = maxiter;
      }
  }
  std::cout << max_elem << std::endl;
  A.diag().print("Eigenvalues A");
  arma::normalise(R).print("Eigenvectors R");

  return 0;
}