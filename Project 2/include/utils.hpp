#ifndef __utils_hpp__  
#define __utils_hpp__

#include <armadillo>
#include <sstream>
#include <string>
#include <iomanip>


// make tridiagonal matrix A = tridag(a, d, a) where A is NxN
arma::mat make_tridiag(int N, double a, double d);

// Return a string with a double in scientific notation
std::string scientific_format(double d, const int width=20, const int prec=10);

// Return a string with an armadillo vector in scientific notation
std::string scientific_format(const arma::vec& v, const int width=20, const int prec=10);

// write to file
void write_file(std::string filename, std::vector< std::vector<double> > values);

//find analytic solution to tridiagonal matrix A = tridiag(a,d,a) where A in NxN
void solve_analytic(int N, double a, double d, arma::vec& eigvals, arma::mat& eigvecs);

//function for comparing eigenvecs+eigenvals with analytic solution
void compare_solutions(arma::mat eigvecs_analytic, arma::vec eigvals_analytic, arma::mat eigvecs, arma::vec eigvals, double eps);

#endif 