#include <iostream>   
#include <string> 
#include <sstream>
#include <cmath>
#include <armadillo>

#include "Particle.hpp"
#include "PenningTrap.hpp"

void simulate_n_particles(
    std::vector<double>&results, 
    PenningTrap penning_trap, 
    int n, double tot_time, 
    double dt, 
    double f, 
    double w_V, 
    std::string evolve_method,
    bool are_interacting
    );
void run_resonance_experiment(
    double f, // amplitude
    bool are_interacting, // if there is Coloumb interaction between particles
    int n, // num. partivles to simulate
    double tot_time, 
    double dt,
    double w_V_start,
    double w_V_stop, 
    double dw_V, 
    std::string evolve_method,
    double B0 = 9.65 * 10,
    double V0 = 2.41 * 1000000,
    double d = 500
    );


int main()
{
    arma::arma_rng::set_seed(2022);

    bool are_interacting = false;
    int n = 100;
    double tot_time = 500;
    double dt = 0.05;

    // frequency
    double w_V_start = 0.2;
    double w_V_stop = 2.5;
    double dw_V = 0.02;

    //amplitude
    double f = 0.1;
    run_resonance_experiment(f, are_interacting, n,tot_time, dt, w_V_start, w_V_stop, dw_V,"RK4");

    f = 0.4;
    run_resonance_experiment(f, are_interacting, n,tot_time, dt, w_V_start, w_V_stop, dw_V,"RK4");

    f = 0.7;
    run_resonance_experiment(f, are_interacting, n,tot_time, dt, w_V_start, w_V_stop, dw_V,"RK4");

    return 0;
}

void simulate_n_particles(
    std::vector<double>& results, 
    PenningTrap penning_trap, 
    int n, 
    double tot_time, 
    double dt, 
    double f, 
    double w_V, 
    std::string evolve_method,
    bool are_interacting
    )
{
    penning_trap.add_particle(n);
    penning_trap.particle_interaction(are_interacting);

    int N_trapped_init = penning_trap.num_particles_in_trap();

    penning_trap.set_time_dependent_V(f, w_V);

    double time = 0;
    int N_steps = tot_time/dt;

    // compute position with Forward Euler
    for(int i = 1; i < N_steps; i++)
    {
        time += dt;
        
        if(evolve_method=="RK4")
        {
            penning_trap.evolve_RK4(dt, time);
        } 
        if(evolve_method=="FE")
        {
            penning_trap.evolve_forward_Euler(dt, time);
        }
    }

    int N_trapped_end = penning_trap.num_particles_in_trap();
    double frac = (double) N_trapped_end / (double) N_trapped_init;
    results.push_back(frac);   
}

void run_resonance_experiment(
    double f, // amplitude
    bool are_interacting, // if there is Coloumb interaction between particles
    int n, // num. partivles to simulate
    double tot_time, 
    double dt,
    double w_V_start,
    double w_V_stop, 
    double dw_V, 
    std::string evolve_method,
    double B0,
    double V0,
    double d
    )
{
    // frequency 
    int N_w_V = (w_V_stop - w_V_start)/dw_V; // num. frequencies to test

    // for saving results given amplitude f_1=0.1
    std::vector<double> results;

    std::ostringstream oss;
    oss << "_interaction_" << are_interacting << "_dw_" << dw_V << "_f_" << f << "_tot_time_" << tot_time;
    std::string experiment_info = oss.str();

    std::string filename = "data/fraction" + experiment_info + ".txt";
    std::ofstream ofile;
    ofile.open(filename);
    int width = 15;
    int prec  = 8;
    
    double w_V = w_V_start;    
    for(int i = 0; i < N_w_V; i++)
    {
        std::cout << "i = " << i << std::endl;
        PenningTrap penning_trap(B0, V0, d);
        simulate_n_particles(results, penning_trap, n, tot_time, dt, f, w_V, evolve_method, are_interacting);

        ofile   << std::setw(width) << std::setprecision(prec) << std::scientific << results.at(i) 
                << std::setw(width) << std::setprecision(prec) << std::scientific << w_V 
                << std::endl;

        w_V += dw_V;
    }

    ofile.close();
    std::cout << "f = " << f << " completed." << std::endl;   
}