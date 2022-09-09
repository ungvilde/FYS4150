# include "algorithms.hpp"

#include <iostream>
#include <vector>
#include <cmath>

std::vector< std::vector<double> > general_algorithm(
    std::vector<double> a, // lower diagonal, size N-1
    std::vector<double> b, // main diagonal, size N
    std::vector<double> c, // upper diagonal, size N-1
    std::vector<double> g, // f(x)*h^2, size N
    int N, // num steps
    double h // stepsize
    )
{
    std::vector<double> v(N); // vector for approx uvalues
    std::vector<double> xvalues(N); // xvalues ranging from h to N*h

    double alpha; // variable used to compute a_i/b_{i-1}
    xvalues[0] = h;
    // execute general algorithm
    for(int i=1; i<N; i++)
    {
        // forward substitution
        alpha = a[i]/b[i-1]; // 1 FLOP
        b[i] = b[i] - alpha*c[i-1]; //2 FLOPs
        g[i] = g[i] - alpha*g[i-1]; //2 FLOPs
        xvalues[i] = xvalues[i-1] + h;
    }

    v[N-1] = g[N-1]/b[N-1]; // 1 FLOPs

    for(int i=N-2; i>=0; i--)
    {
        // backward substitution
        v[i] = (g[i] - c[i]*v[i+1])/b[i]; //3 FLOPs
    }

    //Total: 8*n + 1 FLOPs

    std::vector< std::vector<double> > approxvalues;
    approxvalues.push_back(xvalues);
    approxvalues.push_back(v);

    return approxvalues;
}
