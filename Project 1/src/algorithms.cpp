# include "algorithms.hpp"
# include "utils.hpp"

#include <iostream>
#include <vector>
#include <cmath>

std::vector< std::vector<double> > general_algorithm(
    int N // num steps 
    )
{
    std::vector<double> b(N-1, 2.);
    std::vector<double> a(N-2, -1.);
    std::vector<double> c(N-2, -1.);

    std::vector<double> v(N+1); // vector for approx uvalues
    std::vector<double> xvalues(N+1); // xvalues ranging form 0 to 1
    std::vector<double> g(N+1);

    double h = 1./N; //stepsize
    double alpha; // variable used to compute a_i/b_{i-1}
    
    xvalues[0] = 0.;
    xvalues[N] = 1.;
    v[0] = 0.; //boundary conditions
    v[N] = 0.;
    g[0] = 0.; //boundary conditions
    g[N] = 0.;

    // compute rhs of equation    
    for(int i=1; i<N; i++) 
    {
        xvalues[i] = xvalues[i-1] + h;
        g[i] = h*h*f( xvalues[i] ); 
    }
    
    // execute general algorithm
    for(int i=1; i<N; i++)
    {
        // forward substitution
        alpha = a[i]/b[i-1]; // 1 FLOP
        b[i] = b[i] - alpha*c[i-1]; //2 FLOPs
        g[i] = g[i] - alpha*g[i-1]; //2 FLOPs
    }

    v[N-1] = g[N-1]/b[N-2]; // 1 FLOPs

    for(int i=N-2; i>=1; i--)
    {
        // backward substitution
        v[i] = (g[i] - c[i-1]*v[i+1])/b[i-1]; //3 FLOPs
    }

    //Total: 8*n + 1 FLOPs
    std::vector< std::vector<double> > approxvalues;
    approxvalues.push_back(xvalues);
    approxvalues.push_back(v);

    return approxvalues;
}

std::vector< std::vector<double>> special_algorithm(int N)
{
    std::vector<double> b(N-1, 2.);
    std::vector<double> a(N-2, -1.);
    std::vector<double> c(N-2, -1.);

    std::vector<double> v(N+1); // vector for approx uvalues
    std::vector<double> xvalues(N+1); // xvalues ranging form 0 to 1
    std::vector<double> g(N+1);

    double h = 1./N; //stepsize
    
    xvalues[0] = 0.;
    xvalues[N] = 1.;
    v[0] = 0.; //boundary conditions
    v[N] = 0.;
    g[0] = 0.; //boundary conditions
    g[N] = 0.;

    //initialize:
    // compute rhs of equation    
    for(int i=1; i<N; i++) 
    {
        xvalues[i] = xvalues[i-1] + h;
        g[i] = h*h*f( xvalues[i] ); 
    }

    // compute b values for optimization
    for(int i = 0; i<N; i++)
    {
        b[i] = (i+2.)/(i+1.);
    }

    // execute special algorithm
    for(int i=1; i<N; i++)
    {
        // forward substitution
        g[i] = g[i] + g[i-1]/b[i-1]; //2 FLOPs
    }

    v[N-1] = g[N-1]/b[N-1]; // 1 FLOPs

    for(int i=N-2; i>=1; i--)
    {
        // backward substitution
        v[i] = (g[i] + v[i+1])/b[i-1]; //2 FLOPs
    }

    std::vector< std::vector<double> > approxvalues;
    approxvalues.push_back(xvalues);
    approxvalues.push_back(v);

    return approxvalues;
}