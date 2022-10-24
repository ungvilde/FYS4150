#include <iostream>   
#include <string> 
#include <sstream>
#include <cmath>
#include <armadillo>
#include <time.h>

#include "Particle.hpp"
#include "PenningTrap.hpp"

// run experiment with 100 randomly initiated particles
void simulate_n_particles(std::vector<double>&results, PenningTrap penning_trap, int n, double tot_time, 
    double dt, double f, double w_V, std::string evolve_method, bool are_interacting, bool save_position=false);

// run experiment with varying amplitude and resonance
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
    int k, // iteration of experiment
    double B0 = 9.65 * 10,
    double V0 = 2.41 * 1000000,
    double d = 500
    );

int main()
{
    arma::arma_rng::set_seed(2022);

    int n = 100;
    double tot_time = 500;
    double dt = 0.05;

    // frequencies to explore
    double w_V_start = 0.2;
    double w_V_stop = 2.5;
    double dw_V = 0.02;

    // broad search, test varying amplitude and frequencies
    bool are_interacting = false;
    double f = 0.1;
    run_resonance_experiment(f, are_interacting, n,tot_time, dt, w_V_start, w_V_stop, dw_V, "RK4", 0);

    f = 0.4;
    run_resonance_experiment(f, are_interacting, n,tot_time, dt, w_V_start, w_V_stop, dw_V, "RK4", 0);

    f = 0.7;
    run_resonance_experiment(f, are_interacting, n,tot_time, dt, w_V_start, w_V_stop, dw_V, "RK4", 0);

    // close-grained search of frequencies for comparing with and without Coloumb interaction
    w_V_start = 1.3;
    w_V_stop = 1.5;
    //w_V_start = 2.1;
    //w_V_stop = 2.3;
    dw_V = 0.005;
    f = 0.1;
    for(int k=0; k<5; k++)
    {
    are_interacting = true;
    std::cout << "Simulating 100 particles with interaction." << std::endl;
    clock_t t1 = clock();
    run_resonance_experiment(f, are_interacting, n, tot_time, dt, w_V_start, w_V_stop, dw_V, "RK4", k);
    clock_t t2 = clock();
    double duration_seconds = ((double) (t2 - t1)) / CLOCKS_PER_SEC;
    std::cout << "Simulation over in " << duration_seconds << " seconds." << std::endl;

    are_interacting = false;
    run_resonance_experiment(f, are_interacting, n, tot_time, dt, w_V_start, w_V_stop, dw_V, "RK4", k);
    }
    // // get position of a particle with and without interaction
    // double w_V = 1.5;
    // double f = 0.7;
    // are_interacting = true;
    // PenningTrap penning_trap(9.65 * 10, 2.41 * 1000000, 500);
    // double q = 1;
    // double m = 40.078; // atomic mass of Ca+ ion
    // arma::vec r = arma::vec(3).randn() * 0.1 * 500.;
    // arma::vec v = arma::vec(3).randn() * 0.1 * 500.;
    // Particle p(q, m, r, v);
    // penning_trap.add_particle(p);
    // n = 99;
    // std::vector<double> results;

    // std::cout << "simulating with interaction" << std::endl;
    // clock_t t1 = clock();
    // simulate_n_particles(results, penning_trap, n, tot_time, dt, f, w_V, "RK4", are_interacting, true);
    // clock_t t2 = clock();
    // double duration_seconds = ((double) (t2 - t1)) / CLOCKS_PER_SEC;
    // std::cout << "Simulation over in " << duration_seconds << " seconds." << std::endl;

    // are_interacting = false;
    // std::cout << "simulating without interaction" << std::endl;
    // t1 = clock();
    // simulate_n_particles(results, penning_trap, n, tot_time, dt, f, w_V, "RK4", are_interacting, true);
    // duration_seconds = ((double) (t2 - t1)) / CLOCKS_PER_SEC;
    // std::cout << "Simulation over in " << duration_seconds << " seconds." << std::endl;

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
    bool are_interacting,
    bool save_position
    )
{
    penning_trap.add_particle(n);
    penning_trap.particle_interaction(are_interacting);

    int N_trapped_init = penning_trap.num_particles_in_trap();

    penning_trap.set_time_dependent_V(f, w_V);
    double time = 0;
    int N_steps = tot_time/dt;

    arma::vec x_vals = arma::vec(N_steps, arma::fill::zeros);
    arma::vec y_vals = arma::vec(N_steps, arma::fill::zeros);
    arma::vec z_vals = arma::vec(N_steps, arma::fill::zeros);
    // initial values
    x_vals(0) = penning_trap.p[0].position()(0);
    y_vals(0) = penning_trap.p[0].position()(1);
    z_vals(0) = penning_trap.p[0].position()(2);

    // compute position of single particle over time
    for(int i = 1; i < N_steps; i++)
    {
        time += dt;
        
        if(evolve_method=="RK4")
        {
            penning_trap.evolve_RK4(dt, time);
            arma::vec r = penning_trap.p[0].position();
            x_vals(i) = r(0);
            y_vals(i) = r(1);
            z_vals(i) = r(2);
        } 
        if(evolve_method=="FE")
        {
            penning_trap.evolve_forward_Euler(dt, time);
            arma::vec r = penning_trap.p[0].position();
            x_vals(i) = r(0);
            y_vals(i) = r(1);
            z_vals(i) = r(2);            
        }
    }

    int N_trapped_end = penning_trap.num_particles_in_trap();
    double frac = (double) N_trapped_end / (double) N_trapped_init;
    results.push_back(frac);   

    if(save_position)
    {
        std::ostringstream oss;
        oss << "_interaction_" << are_interacting << "_dt_" << dt <<"_wV_" << w_V << "_f_" << f;
        std::string experiment_info = oss.str();

        x_vals.save("data/x_values_resonance" + experiment_info + ".txt", arma::raw_ascii);
        y_vals.save("data/y_values_resonance" + experiment_info + ".txt", arma::raw_ascii);
        z_vals.save("data/z_values_resonance" + experiment_info + ".txt", arma::raw_ascii);
    }
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
    int k,
    double B0,
    double V0,
    double d
    )
{
    int N_w_V = (w_V_stop - w_V_start)/dw_V; // num. frequencies to test

    // for saving results
    std::vector<double> results;

    std::ostringstream oss;
    oss << "_interaction_" << are_interacting << "_dw_" << dw_V << "_f_" << f << "_tot_time_" << tot_time << "_" << k;
    std::string experiment_info = oss.str();

    std::string filename = "data/fraction_2.1_2.3_" + experiment_info + ".txt";
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
