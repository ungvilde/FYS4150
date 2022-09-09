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
std::vector< std::vector<double> > general_algorithm(
    std::vector<double> a, // lower diagonal, size N-1
    std::vector<double> b, // main diagonal, size N
    std::vector<double> c, // upper diagonal, size N-1
    std::vector<double> g, // RHS: f(x)*h^2, size N
    int N, // num steps
    double h
    );

#endif  // end of include guard __utils_hpp__