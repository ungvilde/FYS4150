#include <iostream> 
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <fstream>

#include "utils.hpp"

double f(double x)
{
    return 100.*exp(-10.*x);
}

int main()
{
    int N=100;
    double h=1./N;
    std::vector<double> a(N);
    std::vector<double> b(N-1);
    std::vector<double> c(N-1);
    std::vector<double> g(N);
    double x = 0.;
    for(int i=0; i<N; i++)
    {
        a[i] = 2;
        b[i] = -1;
        c[i] = -1;
        x += h;
        g[i] = f(x);

    }
    std::vector<double> g(N);
    std::vector<double> v(N);

    g[0] = 0;
    for (int i=0; i<N; i++)
    {
        // forward substitution
        
    }


    return 0;
}