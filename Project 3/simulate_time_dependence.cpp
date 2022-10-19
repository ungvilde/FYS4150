#include <iostream>   
#include <string> 
#include <sstream>
#include <cmath>
#include <armadillo>


#include "Particle.hpp"
#include "PenningTrap.hpp"

void simulate_n_particles(std::vector<double>&results, PenningTrap penning_trap, int n, double tot_time, double dt, double f, double w_V, std::string evolve_method);

int main()
{
    arma::arma_rng::set_seed(2022);

    // // Values for constructing Penning trap
    // double B0 = 9.65 * 10;
    // double V0 = 2.41 * 1000000;
    // double d = 500;

    // // Values for simulation
    // double tot_time = 50; 
    // double dt = 0.05;
    // int n = 100; // num. particles

    // // frequency 
    // double w_V_start = 0.2;
    // double w_V_stop = 2.5;
    // double dw_V = 0.02; // step size along w_V axis
    // int N_w_V = (w_V_stop - w_V_start)/dw_V; // num. frequencies to test

    // // for saving results given amplitude f_1=0.1
    // std::vector<double> results_f1;
    // std::string filename = "data/fraction_f_0.1_500us.txt";
    // std::ofstream ofile;
    // ofile.open(filename);
    // int width = 15;
    // int prec  = 8;
    
    // double f1 = 0.1;
    // double w_V = w_V_start;
    
    // for(int i =0; i<N_w_V; i++)
    // {
    //     std::cout << "i = " << i << std::endl;
    //     PenningTrap penning_trap(B0, V0, d);
    //     simulate_n_particles(results_f1, penning_trap, n, tot_time, dt, f1, w_V, "RK4");
    //     ofile   << std::setw(width) << std::setprecision(prec) << std::scientific << results_f1.at(i) 
    //             << std::setw(width) << std::setprecision(prec) << std::scientific << w_V 
    //             << std::endl;


    //     w_V += dw_V;
    // }
    // ofile.close();

    // std::cout << "f1 = 0.1 completed" << std::endl;
    
    // // for saving results given amplitude 0.4
    // std::vector<double> results_f2;
    // filename = "data/fraction_f_0.4_500us.txt";
    // ofile.open(filename);
    
    // double f2 = 0.4;
    // w_V = w_V_start;

    // for(int i =0; i<N_w_V; i++)
    // {
    //     std::cout << "i = " << i << std::endl;
    //     PenningTrap penning_trap(B0, V0, d);
    //     simulate_n_particles(results_f2, penning_trap, n, tot_time, dt, f2, w_V, "RK4");
    //     ofile   << std::setw(width) << std::setprecision(prec) << std::scientific << results_f2.at(i) 
    //             << std::setw(width) << std::setprecision(prec) << std::scientific << w_V << std::endl;

    //     w_V += dw_V;
    // }
    // ofile.close();

    // std::cout << "f2 = 0.4 completed" << std::endl;

    // // for saving results given amplitude 0.7
    // std::vector<double> results_f3;
    // filename = "data/fraction_f_0.7_500us.txt";
    // ofile.open(filename);
    
    // double f3 = 0.7;
    // w_V = w_V_start;

    // for(int i=0; i<N_w_V; i++)
    // {
    //     std::cout << "i = " << i << std::endl;
    //     PenningTrap penning_trap(B0, V0, d);
    //     simulate_n_particles(results_f3, penning_trap, n, tot_time, dt, f3, w_V, "RK4");
    //     ofile   << std::setw(width) << std::setprecision(prec) << std::scientific << results_f3.at(i) 
    //             << std::setw(width) << std::setprecision(prec) << std::scientific << w_V << std::endl;

    //     w_V += dw_V;
    // }
    // ofile.close();

    // std::cout << "f3 = 0.7 completed" << std::endl;

    // Values for constructing Penning trap
    double B0 = 9.65 * 10;
    double V0 = 2.41 * 1000000;
    double d = 500;

    // Values for simulation
    double tot_time = 50; 
    double dt = 0.05;
    int n = 100; // num. particles

    // frequency 
    double w_V_start = 1.0;
    double w_V_stop = 1.7;
    double dw_V = 0.005; // step size along w_V axis
    int N_w_V = (w_V_stop - w_V_start)/dw_V; // num. frequencies to test

    // for saving results given amplitude f_1=0.1
    std::vector<double> results_f1;
    std::string filename = "data/fraction_f_0.1_500us_dw_0.005.txt";
    std::ofstream ofile;
    ofile.open(filename);
    int width = 15;
    int prec  = 8;
    
    double f1 = 0.1;
    double w_V = w_V_start;
    
    for(int i =0; i<N_w_V; i++)
    {
        std::cout << "i = " << i << std::endl;
        PenningTrap penning_trap(B0, V0, d);
        simulate_n_particles(results_f1, penning_trap, n, tot_time, dt, f1, w_V, "RK4");
        ofile   << std::setw(width) << std::setprecision(prec) << std::scientific << results_f1.at(i) 
                << std::setw(width) << std::setprecision(prec) << std::scientific << w_V 
                << std::endl;


        w_V += dw_V;
    }
    ofile.close();

    std::cout << "f1 = 0.1 completed" << std::endl;
    
    // for saving results given amplitude 0.4
    std::vector<double> results_f2;
    filename = "data/fraction_f_0.4_500us_dw_0.005.txt";
    ofile.open(filename);
    
    double f2 = 0.4;
    w_V = w_V_start;

    for(int i =0; i<N_w_V; i++)
    {
        std::cout << "i = " << i << std::endl;
        PenningTrap penning_trap(B0, V0, d);
        simulate_n_particles(results_f2, penning_trap, n, tot_time, dt, f2, w_V, "RK4");
        ofile   << std::setw(width) << std::setprecision(prec) << std::scientific << results_f2.at(i) 
                << std::setw(width) << std::setprecision(prec) << std::scientific << w_V << std::endl;

        w_V += dw_V;
    }
    ofile.close();

    std::cout << "f2 = 0.4 completed" << std::endl;

    // for saving results given amplitude 0.7
    std::vector<double> results_f3;
    filename = "data/fraction_f_0.7_500us_dw_0.005.txt";
    ofile.open(filename);
    
    double f3 = 0.7;
    w_V = w_V_start;

    for(int i=0; i<N_w_V; i++)
    {
        std::cout << "i = " << i << std::endl;
        PenningTrap penning_trap(B0, V0, d);
        simulate_n_particles(results_f3, penning_trap, n, tot_time, dt, f3, w_V, "RK4");
        ofile   << std::setw(width) << std::setprecision(prec) << std::scientific << results_f3.at(i) 
                << std::setw(width) << std::setprecision(prec) << std::scientific << w_V << std::endl;

        w_V += dw_V;
    }
    ofile.close();

    std::cout << "f3 = 0.7 completed" << std::endl;

    return 0;
}

void simulate_n_particles(std::vector<double>& results, PenningTrap penning_trap, int n, double tot_time, double dt, double f, double w_V, std::string evolve_method)
{
    penning_trap.add_particle(n);
    penning_trap.particle_interaction(false);

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
