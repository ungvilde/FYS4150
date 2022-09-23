#include <iostream>
#include <cmath>
#include "utils.hpp"
#include "algorithm.hpp"

int main()
{
  // Problem 2: 
  // Compare armadillo and analyitcal solution
  int N = 6; // number of equations to solve
  double h = 0.01; //stepsize
  double d = 2. / (h*h); // main diagonal
  double a = -1. / (h*h); // upper + lower diagonal

  // make tridiagonal matrix
  arma::mat A = make_tridiag(N, a, d);

  // solve eigenvalue problem using armadillo
  arma::vec eigvals_arma;
  arma::mat eigvecs_arma;
  arma::eig_sym(eigvals_arma, eigvecs_arma, A);
  eigvecs_arma = arma::normalise(eigvecs_arma); // normalise for comparison with analytic solution
  
  // find analytic solution
  arma::vec eigvals_analytic(N); //eigenvalues
  arma::mat eigvecs_analytic = arma::mat(N,N); //eigenvectors
  solve_analytic(N, a, d, eigvals_analytic, eigvecs_analytic);

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
  A_test.print("A_test:"); 
  
  int l;
  int k;
  double max_elem = max_offdiag_symmetric(A_test, k, l);
  std::cout << "Max. absolute off-diagonal element: " << max_elem << std::endl;
  std::cout << "With indices: k=" << k << " and l=" << l << std::endl;

  // Problem 4: 
  // test jacobi rotation algorithm
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  int iterations;
  int maxiter = 100;
  bool converged = false;

  A = make_tridiag(6, a, d);
  jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);
  std::cout << "Converged after " << iterations << " iterations." << std::endl;

  eigenvalues.print("Jacobi rotation eigenvalues: ");
  eigenvectors.print("Jacobi rotation eigenvectors: ");

  std::cout << "Comparing results with analytical solution:" << std::endl;
  compare_solutions(eigvecs_analytic, eigvals_analytic, eigenvectors, eigenvalues, eps);

  // Problem 5:
  // checking how num iterations of algorithm scales with matrix size
  // int M = 8; // num. matrices of varying N to compute
  // arma::mat data_values = arma::mat(M, 3);
  // N = 4;
  // // loop through different possibilities for N
  // for(int i=0; i < M; i++){
  //   std::cout << "N = " << N << std::endl; 
  //   A = make_tridiag(N, a, d);
  //   maxiter = 1000000;
  //   iterations = 0;
  //   converged = false;
  //   arma::vec eigenvalues_N;
  //   arma::mat eigenvectors_N;
  //   jacobi_eigensolver(A, eps, eigenvalues_N, eigenvectors_N, maxiter, iterations, converged);
  //   data_values(i, 0) = N;
  //   data_values(i, 1) = iterations;
  //   data_values(i, 2) = converged;
  //   std::cout << "For N = " << N << std::endl;
  //   std::cout << "Num. iterations: " << iterations << " Converged? " << converged << " (Yes=1/No=0)" << std::endl;
  //   N = 2*N; // update num. equations to solve
  // }

  // // find num iterations (and if it converged) as function of N
  // data_values.save("data/problem5.txt", arma::raw_ascii);
  // from here, make plot in python with plot.py

  // problem 6
  // find results for N=10
  N = 10;
  A = make_tridiag(N, a, d);

  jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);
  std::cout << "Converged after " << iterations << " iterations." << std::endl;

  eigenvalues.print("Jacobi rotation eigenvalues: ");
  eigenvectors.print("Jacobi rotation eigenvectors: ");  
  arma::uvec inds = { 0, 1, 2 }; // get the three eigenvectors corresponding to the three lowest eigenvalues
  eigenvectors.cols(inds).print("Solutions"); //save("data/problem6.txt", arma::raw_ascii)
  arma::mat V = eigenvectors.cols(inds);
  V.save("data/problem6_numerical.txt", arma::raw_ascii);

  arma::vec eigvals_analytic6(N); //eigenvalues
  arma::mat eigvecs_analytic6 = arma::mat(N,N); //eigenvectors
  solve_analytic(N, a, d, eigvals_analytic6, eigvecs_analytic6);
  arma::mat V_analytic = eigvecs_analytic6.cols(inds);
  V_analytic.save("data/problem6_analytic.txt", arma::raw_ascii);
  return 0;
}