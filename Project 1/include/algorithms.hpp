//
// A collection of "clever" algorithms
//

// First we add an "include guard". It ensures that this header file can 
// only end up being included *once* for the compilation of any given .cpp file
#ifndef __algorithms_hpp__  
#define __algorithms_hpp__

#include <iostream>
#include <vector>

// algorithm for solving general tridiagonal matrix equation
std::vector< std::vector<double> > general_algorithm(int N);

// algorithm for solving special case of tridiagonal matrix eq. where a
std::vector< std::vector<double>> special_algorithm(int N);

// script for timing respective algorithms
double time_general_algorithm(int N, int num_iter);

double time_special_algorithm(int N, int num_iter);


#endif  // end of include guard __utils_hpp__