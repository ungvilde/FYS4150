#include <iostream>
#include <cmath>
#include "utils.hpp"
#include "algorithm.hpp"

int main()
{
  // Problem 2: 
  // Compare armadillo and analytical solution
  int N = 6; // number of equations to solve
  double h = 1/N; //stepsize
  double d = 2. / (h*h); // main diagonal
  double a = -1. / (h*h); // upper + lower diagonal

  // make tridiagonal matrix
  arma::mat A = make_tridiag(N);

  // solve eigenvalue problem using armadillo
  arma::vec eigvals_arma;
  arma::mat eigvecs_arma;
  arma::eig_sym(eigvals_arma, eigvecs_arma, A);
  eigvecs_arma = arma::normalise(eigvecs_arma); // normalise for comparison with analytic solution
  
  // find analytic solution
  arma::vec eigvals_analytic(N); //eigenvalues
  arma::mat eigvecs_analytic = arma::mat(N,N); //eigenvectors
  solve_analytic(N, eigvals_analytic, eigvecs_analytic);

  // now we compare the solutions
  double eps = 1e-10;
  std::cout << "Solution to problem 2:" << std::endl;
  compare_solutions(eigvecs_analytic, eigvals_analytic, eigvecs_arma, eigvals_arma, eps);

  // Problem 3: 
  // test max_offdiag_symmetric 
  std::cout << "Solution to problem 3:" << std::endl;
  arma::mat A_test = arma::mat(4,4).fill(0.);
  A_test.diag() = arma::vec(4).fill(1.);
  A_test(3,0) = 0.5;
  A_test(2,1) = -0.7;
  A_test += A_test.t();
  A_test.print("Matrix for testing:"); 
  
  int l;
  int k;
  double max_elem = max_offdiag_symmetric(A_test, k, l);
  std::cout << "Max. absolute off-diagonal element: " << max_elem << std::endl;
  std::cout << "With indices: k=" << k << " and l=" << l << std::endl;

  // Problem 4: 
  // test jacobi rotation algorithm
  arma::vec eigenvals_numerical;
  arma::mat eigenvecs_numerical;
  int iterations = 0;
  int maxiter = 1000;
  bool converged = false;

  A = make_tridiag(6);
  A.print("A matrix");
  jacobi_eigensolver(A, eps, eigenvals_numerical, eigenvecs_numerical, maxiter, iterations, converged);
  std::cout << "Converged after " << iterations << " iterations." << std::endl;

  eigenvals_numerical.print("Jacobi rotation eigenvalues: ");
  eigenvecs_numerical.print("Jacobi rotation eigenvectors: ");

  std::cout << "Comparing results with analytical solution:" << std::endl;
  compare_solutions(eigvecs_analytic, eigvals_analytic, eigenvecs_numerical, eigenvals_numerical, eps);

  // Problem 5:
  // checking how num. iterations of algorithm scales with matrix size
  int M = 7; // num. solutions of varying N to compute
  arma::mat data_values = arma::mat(M, 3);
  N = 4; // start with 2^2 (up to 2^8)
  maxiter = 1000000;

  // loop through different possibilities for N
  for(int i=0; i < M; i++){
    std::cout << "N = " << N << std::endl; 
    A = make_tridiag(N);
    iterations = 0;
    converged = false;
    arma::vec eigenvalues_N;
    arma::mat eigenvectors_N;

    jacobi_eigensolver(A, eps, eigenvalues_N, eigenvectors_N, maxiter, iterations, converged);
    data_values(i, 0) = N;
    data_values(i, 1) = iterations;
    data_values(i, 2) = converged;

    std::cout << "For N = " << N << std::endl;
    std::cout << "Num. iterations: " << iterations << " Converged? " << converged << " (Yes=1/No=0)" << std::endl;
    N = 2*N; // update num. equations to solve
  }

  data_values.save("data/problem5.txt", arma::raw_ascii);
  // from here, make plot in python with plot.py

  // problem 6
  // find and plot results for N=10
  N = 10;
  A = make_tridiag(N);
  arma::vec eigenvalues_N10;
  arma::mat eigenvectors_N10;

  jacobi_eigensolver(A, eps, eigenvalues_N10, eigenvectors_N10, maxiter, iterations, converged);
  std::cout << "Converged after " << iterations << " iterations." << std::endl;

  arma::uvec inds = { 0, 1, 2 }; // get the three eigenvectors corresponding to the three lowest eigenvalues
  eigenvectors_N10.cols(inds).print("Solutions"); 
  arma::mat V = eigenvectors_N10.cols(inds);
  arma::mat lambda = eigenvalues_N10(inds);
  V.save("data/problem6_numerical_N10.txt", arma::raw_ascii);
  lambda.save("data/problem6_numerical_N10_eigenvals.txt", arma::raw_ascii);

  arma::vec eigvals_analytic_N10(N); //eigenvalues
  arma::mat eigvecs_analytic_N10 = arma::mat(N,N); //eigenvectors
  solve_analytic(N, eigvals_analytic_N10, eigvecs_analytic_N10);
  arma::mat V_analytic_N10 = eigvecs_analytic_N10.cols(inds);
  V_analytic_N10.save("data/problem6_analytic_N10.txt", arma::raw_ascii);

  // repeat for N=100
  N = 100;
  A = make_tridiag(N);
  arma::vec eigenvalues_N100;
  arma::mat eigenvectors_N100;

  jacobi_eigensolver(A, eps, eigenvalues_N100, eigenvectors_N100, maxiter, iterations, converged);
  std::cout << "Converged after " << iterations << " iterations." << std::endl;
  V = eigenvectors_N100.cols(inds);
  lambda = eigenvalues_N100(inds);
  V.save("data/problem6_numerical_N100.txt", arma::raw_ascii);
  lambda.save("data/problem6_numerical_N100_eigenvals.txt", arma::raw_ascii);

  arma::vec eigvals_analytic_N100(N); //eigenvalues
  arma::mat eigvecs_analytic_N100 = arma::mat(N,N); //eigenvectors
  solve_analytic(N, eigvals_analytic_N100, eigvecs_analytic_N100);
  arma::mat V_analytic_N100 = eigvecs_analytic_N100.cols(inds);
  V_analytic_N100.save("data/problem6_analytic_N100.txt", arma::raw_ascii);

  return 0;
}