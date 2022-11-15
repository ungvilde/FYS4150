#include <iostream>   
#include <string> 
#include <sstream>
#include <cmath>

#include "IsingModel.hpp"


// CPATH=/opt/homebrew/include LIBRARY_PATH=/opt/homebrew/lib g++ main.cpp -std=c++11 -larmadillo src/* -I include -o main

int main()
{
    int L = 2;

    IsingModel ising_model(2, 1.0);
    ising_model.run_n_MC_cycles(50000, 1000);

    double Z = 2 * exp(8.0) + 12 + 2*exp(-8.0);
    double analytic_energy = - 1.0 / Z * (2*8 * exp(8.0) - 2*8 * exp(-8.0));
    double analytic_energy2 = 128.0 / Z * (exp(8.0) + exp(-8.0));
    double analytic_heat_capacity = 1.0 / (L*L) * (analytic_energy2 - analytic_energy * analytic_energy);

    std::cout << "Analytic expected energy = " << analytic_energy << std::endl;
    std::cout << "Analytic expected squared energy = " << analytic_energy2 << std::endl;

    std::cout << "Analytic expected heat capacity = " << analytic_heat_capacity << std::endl;

    return 0;
}