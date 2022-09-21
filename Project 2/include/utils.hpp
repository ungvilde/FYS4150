#ifndef __utils_hpp__  
#define __utils_hpp__

#include <armadillo>
#include <sstream>
#include <string>
#include <iomanip>

// Return a string with a double in scientific notation
std::string scientific_format(double d, const int width=20, const int prec=10);

// Return a string with an armadillo vector in scientific notation
std::string scientific_format(const arma::vec& v, const int width=20, const int prec=10);

#endif 