#ifndef __IsingModel_hpp__  
#define __IsingModel_hpp__

#include <random>
#include <chrono>
#include <iostream>
#include <armadillo>
#include <cmath>
#include <stdexcept>


class IsingModel
{

    public:

        // constructor
        IsingModel(int L, double T); 

        // lattice variables
        int L; // lattice size
        arma::mat lattice; // lattice as armadillo matrix object
        double J = 1.0; // coupling
        double kb = 1.0; //boltzman constant
        double T; //temperature, in units of J/kb
        double beta; // 1 / T*kb
        int N; // L*L

        // quantities
        double energy; // energy of lattice
        double magnetisation; // magnetisation

        // method for initialising the lattice
        void initialize_lattice(std::string initialisation); 

        // method for running through a single monte carly cycle
        void run_MC_cycle();

        // method for running through n monte carlo cycles, storing the E, M, etc for each cycle
        // n_cycles: total number of MC cycles
        // n0: determines the burn-in, ie how many iterations we run before we store the data values
        // initialisation: deterimines the configuration of the lattice and the beginning of the simulation
        arma::mat run_n_MC_cycles(int n_cycles, int n0, std::string initialisation); 

        // method for computing the energy of a lattice
        double compute_energy();

        // method for computing magnetisation of a lattice
        double compute_magnetisation();
};

#endif