#ifndef __utils_hpp__  
#define __utils_hpp__

#include <armadillo>
#include <sstream>
#include <string>
#include <iomanip>


// make tridiagonal matrix A = tridag(a, d, a) where A is NxN
arma::mat make_tridiag(int N, double a, double d);

//find analytic solution to tridiagonal matrix A = tridiag(a,d,a) where A in NxN
void solve_analytic(int N, double a, double d, arma::vec& eigvals, arma::mat& eigvecs);

//function for comparing eigenvecs+eigenvals with analytic solution
void compare_solutions(arma::mat eigvecs_analytic, arma::vec eigvals_analytic, arma::mat eigvecs, arma::vec eigvals, double eps);

#endif 