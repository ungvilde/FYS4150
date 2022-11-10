#include <random>
#include <chrono>
#include <iostream>
#include <armadillo>
#include <cmath>

// CPATH=/opt/homebrew/include LIBRARY_PATH=/opt/homebrew/lib g++ main.cpp -std=c++11 -larmadillo -o main

double compute_magnetisation(arma::mat lattice);
double compute_energy(arma::mat lattice);
double compute_heat_capacity(arma::mat lattice);

int main()
{
    double J = 1.0; // coupling constant
    double T = 1.0; // temperature
    double kb = 1.0; // boltzmann constant
    double beta = 1.0 / (kb * T);
    int L = 5; // lattice length
    int N = L*L; // num. spins

    // initialize lattice
    // can be a function (or we could go crazy and make a class, with initialization as a method)
    arma::mat lattice = arma::randi<arma::mat>(L, L, arma::distr_param(0, 1));
    lattice = 2.0 * lattice - 1;
    std::cout << "L = " << L << std::endl;
    lattice.print("Initial lattice:");

    // change in energy of the system, updated at each attempted spin flip
    double dE = 0.0;

    for(int k=0; k<N; k++)
    {
        int i = arma::randi(arma::distr_param(0, L-1));
        int j = arma::randi(arma::distr_param(0, L-1));
        double energy_noflip = - J * lattice(i, j) * (
            lattice(i, ((j - 1) % L + L) % L) +
            lattice(i, ((j + 1) % L + L) % L) +
            lattice(((i - 1) % L + L) % L, j) +
            lattice(((i + 1) % L + L) % L, j)
        );
        
        double energy_flip = -1 * energy_noflip;
        double u = arma::randu(arma::distr_param(0,1));
                
        if(energy_flip - energy_noflip <= 0)
        {
            lattice(i, j) *= -1; // do the flip
            dE += energy_flip;
        }
        else if(u <= exp(beta * (energy_flip - energy_noflip) ))
        {
            lattice(i, j) *= -1; // do the flip
            dE += energy_flip; // update energy change in system
        }
    }

    lattice.print("Final lattice");

    return 0;
}

double compute_magnetisation(arma::mat lattice);

double compute_energy(arma::mat lattice, double J)
{
    int L = lattice.n_cols;
    double energy = 0;
    for(int i = 0; i<L; i++)
    {
        for(int j = 0; i<L; i++)
        {
            energy += -J * lattice(i, j)*(
                lattice(i, ((j - 1) % L + L) % L) +
                lattice(((i + 1) % L + L) % L, j) 
            );
        }

    }
}