// First we add an "include guard". It ensures that this header file can 
// only end up being included *once* for the compilation of any given .cpp file
#ifndef __utils_hpp__  
#define __utils_hpp__

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

// Return a string with a double in scientific notation
std::string scientific_format(double d, const int width=20, const int prec=10);

// write file with uvalues anf xvalues computed at 100 points 
void write_file(std::string filename, std::vector< std::vector<double> > values);

// exact solution
double u(double x);
double f(double x);

#endif  // end of include guard __utils_hpp__

