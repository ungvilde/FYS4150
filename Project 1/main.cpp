#include <iostream> 
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>

#include "utils.hpp"


int main()
{
    int N = 100;
    double h = 1./N;
    std::vector<double> xvalues(N+1);
    std::vector<double> uvalues(N+1);
    xvalues[0] = 0.;
    uvalues[0] = 0.;
    //then loop to assign values incremented by small step size to every element in vector
    for (int i = 1; i < N+1; i++)
    {
        xvalues[i] = xvalues[i-1] + h;
        uvalues[i] = u(xvalues[i]);
    }
    // save in data file
    std::vector< std::vector<double> > datavalues;
    datavalues.push_back(xvalues);
    datavalues.push_back(uvalues);
    write_file("data.txt", datavalues);

    // now execute general algorithm

    return 0;
}

