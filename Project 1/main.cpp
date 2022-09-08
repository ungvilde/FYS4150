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

    // prepare matrix equation
    std::vector<double> g(N);
    //std::vector<double> gtilde(N);

    std::vector<double> b(N, 2.);
    //std::vector<double> btilde(N);

    std::vector<double> a(N-1, -1.);
    std::vector<double> c(N-1, -1.);
    std::vector<double> v(N); // variable to be solved

    for(int i=0; i<N; i++)
    {
        g[i] = h*h*f(xvalues[i]); // compute rhs of equation
    }
    
    double alpha; // variable used to compute c_i/a_{i-1}

    // execute general algorithm
    for(int i=1; i<N; i++)
    {
        // forward substitution
        alpha = a[i]/b[i-1];
        b[i] = b[i] - alpha*c[i-1];
        g[i] = g[i] - alpha*g[i-1];
    }

    v[N-1] = g[N-1]/b[N-1];
    for(int i=N-2; i>=0; i--)
    {
        v[i] = (g[i] - c[i]*v[i+1])/b[i];
    }
    std::vector<std::vector<double>> approxvalues;
    approxvalues.push_back(xvalues);
    approxvalues.push_back(v);
    write_file("approx.txt", approxvalues);

    return 0;
}

