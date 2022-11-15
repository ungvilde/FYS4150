#include <iostream>   
#include <string> 
#include <sstream>
#include <cmath>

#include "IsingModel.hpp"


// CPATH=/opt/homebrew/include LIBRARY_PATH=/opt/homebrew/lib g++ main.cpp -std=c++11 -larmadillo src/* -I include -o main

int main()
{
    int L = 20;

    IsingModel ising_model0(L, 1.0, "random");
    arma::mat data0 = ising_model0.run_n_MC_cycles(50000, 0);
    data0.save("datasets/data_L20_initRand_MC50000_T1.txt", arma::raw_ascii);
    
    IsingModel ising_model1(L, 2.4, "random");
    arma::mat data1 = ising_model1.run_n_MC_cycles(50000, 0);
    data1.save("datasets/data_L20_initRand_MC50000_T2.4.txt", arma::raw_ascii);

    IsingModel ising_model2(L, 1.0, "ordered");
    arma::mat data2 = ising_model2.run_n_MC_cycles(50000, 0);
    data2.save("datasets/data_L20_initOrdered_MC50000_T1.txt", arma::raw_ascii);

    IsingModel ising_model3(L, 2.4, "ordered");
    arma::mat data3 = ising_model3.run_n_MC_cycles(50000, 0);
    data3.save("datasets/data_L20_initOrdered_MC50000_T2.4.txt", arma::raw_ascii);

    // double Z = 2 * exp(8.0) + 12 + 2*exp(-8.0);
    // double analytic_energy = - 1.0 / Z * (2*8 * exp(8.0) - 2*8 * exp(-8.0));
    // double analytic_energy2 = 128.0 / Z * (exp(8.0) + exp(-8.0));
    // double analytic_heat_capacity = 1.0 / (L*L) * (analytic_energy2 - analytic_energy * analytic_energy);

    // std::cout << "Analytic expected energy = " << analytic_energy << std::endl;
    // std::cout << "Analytic expected squared energy = " << analytic_energy2 << std::endl;

    // std::cout << "Analytic expected heat capacity = " << analytic_heat_capacity << std::endl;

    return 0;
}