#ifndef __IsingModel_hpp__  
#define __IsingModel_hpp__

#include <random>
#include <chrono>
#include <iostream>
#include <armadillo>
#include <cmath>


class IsingModel
{

    public:

        // constructor
        IsingModel(int L, double T, std::string initialisation); 

        // lattice variables
        int L; // lattice size
        arma::mat lattice; // lattice as armadillo matrix object
        double J = 1.0; // coupling
        double kb = 1.0; //boltzman constant
        double T; //temperature
        double beta;
        int N; // L*L
        std::string initialisation;

        // measurables
        double energy; // energy of lattice
        double energy2; // squared energy
        double magnetisation; // magnetisation
        double magnetisation_abs; // magnetisation in absolute value

        // // MC parameters
        // int n_cycles; // num. MC cycles
        // int n0; // burntime, ie num. MC cycles to run before we store the values

        

        // for initialising the lattice
        void initialize_lattice(std::string initialisation); 

        // for running through a single monte carly cycle
        void run_MC_cycle();

        // for running through n monte carlo cycles, storing the E, M, etc for each cycle
        void run_n_MC_cycles(int n_cycles, int n0); 

        // for computing the energy of a lattice
        double compute_energy();

        // for computing magnetisation of a lattice
        double compute_magnetisation();
};

#endif