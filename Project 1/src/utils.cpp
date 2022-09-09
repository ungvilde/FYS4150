//
// Defintions for functions declared in utils.hpp
//

// By including utils.hpp, we also include all 
// the headers included in utils.hpp (<armadillo>, <sstream>, ...)

#include "utils.hpp"


void write_file(std::string filename, std::vector< std::vector<double> > values)
{
    // this is a more general and usable version
    // first we create a file and write the data points to this file
    std::ofstream ofile;
    ofile.open(filename);
    
    // width and precision parameters used to format the output
    int width = 20;
    int precision = 10;
    // Loop over steps
    int M = values[0].size();
    for(int j = 0; j < M; j++)
    {
    ofile << scientific_format(values[0][j], width, precision) //xvalues
        << scientific_format(values[1][j], width, precision) //yvalues
        << std::endl;
    }

    // close the data file
    ofile.close();

}

// define exact solution
double u(double x)
{
    return 1 - (1 - exp(-10.))*x - exp(-10*x);
}

double f(double x)
{
    return 100.*exp(-10*x);
}

// Return a string with a double in scientific notation (NOTE: written by Anders)
std::string scientific_format(const double d, const int width, const int prec)
{
  std::stringstream ss;
  ss << std::setw(width) << std::setprecision(prec) << std::scientific << d;
  return ss.str();
}