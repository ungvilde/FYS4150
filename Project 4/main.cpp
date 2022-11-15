#include <iostream>   
#include <string> 
#include <sstream>
#include <cmath>

#include "IsingModel.hpp"

int L; // variable used for lattice geometry
double T;
// CPATH=/opt/homebrew/include LIBRARY_PATH=/opt/homebrew/lib g++ main.cpp -std=c++11 -larmadillo src/* -I include -o main

int main()
{
    // ---------
    // Problem 4
    // ---------
    // L = 2;
    // T = 1.0;
    // double Z = 2 * exp(8.0) + 12 + 2*exp(-8.0); // partition function
    // double analytic_energy = - 1.0 / Z * (16 * exp(8.0) - 16 * exp(-8.0));
    // double analytic_eps = analytic_energy / (L * L);

    // double analytic_energy2 = 128.0 / Z * (exp(8.0) + exp(-8.0));
    // double analytic_eps2 = analytic_energy2 / (L*L * L*L);

    // double analytic_magnetisation = 8.0 / Z * (exp(8.0) + 2);
    // double analytic_magnetisation2 = 32.0 / Z * (exp(8.0) + 1);

    // double analytic_heat_capacity = 1.0 / (L*L) * (analytic_energy2 - analytic_energy * analytic_energy);
    // double analytic_susceptibility = 1.0 / L * (analytic_magnetisation2 - analytic_magnetisation*analytic_magnetisation);

    // std::cout << "Analytic expected energy = " << analytic_energy / (L*L) << std::endl;
    // //std::cout << "Analytic expected squared energy = " << analytic_energy2 << std::endl;
    // std::cout << "Analytic expected heat capacity = " << analytic_heat_capacity << std::endl;
    // std::cout << "Analytic expected abs. magnetisation = " << analytic_magnetisation / (L*L) << std::endl;
    // //std::cout << "Analytic expected squared magnetisation = " << analytic_magnetisation2 << std::endl;
    // std::cout << "Analytic expected susceptibility = " << analytic_susceptibility << std::endl;

    // IsingModel ising_model(L, T);
    // arma::mat data = ising_model.run_n_MC_cycles(10000, 0, "random");
    // std::cout << "Numeric expected energy = " << mean(data.col(0))/(L*L) << std::endl;
    // //std::cout << "Numeric expected squared energy = " << mean( data.col(0) % data.col(0)) << std::endl;
    // std::cout << "Numeric expected heat capacity = " << 1.0/ (L*L) * (mean( data.col(0) % data.col(0)) - mean(data.col(0))*mean(data.col(0))) << std::endl;
    // std::cout << "Numeric expected abs. magnetisation = " << abs(mean(data.col(1))/(L*L)) << std::endl;
    // std::cout << "Numeric expected susceptibility = " << 1.0 / L * (mean( data.col(1) % data.col(1)) - abs(mean(data.col(1)))*abs(mean(data.col(1)))) << std::endl;

    // ----------
    // Problem 5:
    // ----------
    // L = 20;
    // IsingModel ising_model_T1(L, 1.0);
    // arma::mat data0 = ising_model_T1.run_n_MC_cycles(50000, 0, "random");
    // data0.save("datasets/data_L20_initRand_MC50000_T1.txt", arma::raw_ascii);
    // arma::mat data2 = ising_model_T1.run_n_MC_cycles(50000, 0, "ordered");
    // data2.save("datasets/data_L20_initOrdered_MC50000_T1.txt", arma::raw_ascii);

    // IsingModel ising_model_T24(L, 2.4);
    // arma::mat data1 = ising_model_T24.run_n_MC_cycles(50000, 0, "random");
    // data1.save("datasets/data_L20_initRand_MC50000_T2.4.txt", arma::raw_ascii);
    // arma::mat data3 = ising_model_T24.run_n_MC_cycles(50000, 0, "ordered");
    // data3.save("datasets/data_L20_initOrdered_MC50000_T2.4.txt", arma::raw_ascii);

    // ----------
    // Problem 6:
    // ----------

    

    return 0;
}