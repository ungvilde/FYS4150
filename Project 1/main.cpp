// Write a program that defines a vector of values, 
// evaluates the exact solution  above for these points, and 
// outputs the  and  values as two columns in a data file. 
// The numbers should be written to file in scientific notation
// and with a fixed number of decimals. (Choose a sensible number.)

#include <iostream> 
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>


double u(double x);

int main()
{
    int N = 100; //num points to evaluate u at
    double dx = 1./N;
    // create a vector of length N+1
    std::vector<double> xvalues(N+1);
    // assign initial value to the vector directly
    xvalues[0] = 0.;

    //then loop to assign values incremented by small step size to every element in vector
    for (int i = 1; i < N+1; i++)
    {
        xvalues[i] = xvalues[i-1] + dx;
    }

    // compute exact values of function u 
    std::vector<double> uvalues(N+1);
    for (int i = 0; i < N+1; i++)
    {
        uvalues[i] = u(xvalues[i]);
    }

    // now we create a file and write the data points to this file
    std::string filename = "data.txt";

    // create and open the file, and connect it to the filename given
    std::ofstream ofile;
    ofile.open(filename);
    
    // width and precision parameters used to format the output
    int width = 12;
    int precision = 4;

    // Loop over steps
    for (int i = 0; i < N+1; i++)
    {
    // Write a line with the current x and y values (nicely formatted) to file
    ofile << std::setw(width) << std::setprecision(precision) << std::scientific << xvalues[i]
          << std::setw(width) << std::setprecision(precision) << std::scientific << uvalues[i]
          << std::endl;
    }  

    // close the data file
    ofile.close();

    return 0;
}

// define exact solution
double u(double x){
    return 1 - (1 - exp(-10.)*x - exp(-10*x));
}