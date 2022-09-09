#include <iostream> 
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>

#include "utils.hpp"
#include "algorithms.hpp"

int main()
{
    int N = 100; //number of steps along x-axis to compute exact values
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


    // now we make approximation
    int Nsteps = 1000; // number of steps excluding boundary points x=0 and x=1
    h = 1./Nsteps;
    // prepare matrix equation
    std::vector<double> g(Nsteps);
    std::vector<double> b(Nsteps, 2.);
    std::vector<double> a(Nsteps-1, -1.);
    std::vector<double> c(Nsteps-1, -1.);
    std::vector<double> v(Nsteps);
    double x = 0;

    for(int i=0; i<Nsteps; i++)
    {
        g[i] = h*h*f(x); // compute rhs of equation
        x += h;
    }

    std::vector< std::vector<double> > approxvalues = general_algorithm(a,b,c,g, Nsteps, h);  
    approxvalues[0].insert(approxvalues[0].begin(), 0.); // add x=0 and v=0 to include boundry conditions
    approxvalues[1].insert(approxvalues[1].begin(), 0.);
    write_file("approx.txt", approxvalues);

    return 0;
}

