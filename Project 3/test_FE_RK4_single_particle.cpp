#include <iostream>    

#include "Particle.hpp"
#include "PenningTrap.hpp"


int main()
{
    double q = 1;
    double m = 40.078; // atomic mass of Ca+ ion
    arma::vec r1 = arma::vec("20. 0. 20.");
    arma::vec v1 = arma::vec("0. 25. 0.");
    Particle p1(q, m, r1, v1);
    
    double B0 = 9.65 * 10;
    double V0 = 2.41*1000000;
    double d = 500;
    PenningTrap test_trap(B0, V0, d);

    test_trap.add_particle(p1);

    double dt = 0.001; 
    int N_steps = 50000; // 50us/dt

    arma::vec x_vals = arma::vec(N_steps, arma::fill::zeros);
    arma::vec y_vals = arma::vec(N_steps, arma::fill::zeros);
    arma::vec z_vals = arma::vec(N_steps, arma::fill::zeros);
    
    // initial values
    x_vals(0) = 20.;
    y_vals(0) = 0.;
    z_vals(0) = 20.;

    // compute position with Forward Euler
    for(int i = 1; i < N_steps; i++)
    {
        test_trap.evolve_forward_Euler(dt);
        arma::vec r = test_trap.p[0].position();
        x_vals(i) = r(0);
        y_vals(i) = r(1);
        z_vals(i) = r(2);
    }

    x_vals.save("data/x_values_FE.txt", arma::raw_ascii);
    y_vals.save("data/y_values_FE.txt", arma::raw_ascii);
    z_vals.save("data/z_values_FE.txt", arma::raw_ascii);

    PenningTrap test_trap2(B0, V0, d);
    test_trap2.add_particle(p1);
    //test_trap2.add_particle(p2);

    // compute position with Runge-Kutta 4th order
    for(int i = 1; i < N_steps; i++)
    {
        test_trap2.evolve_RK4(dt);
        arma::vec r = test_trap2.p[0].position();
        x_vals(i) = r(0);
        y_vals(i) = r(1);
        z_vals(i) = r(2);
    }

    x_vals.save("data/x_values_RK4.txt", arma::raw_ascii);
    y_vals.save("data/y_values_RK4.txt", arma::raw_ascii);
    z_vals.save("data/z_values_RK4.txt", arma::raw_ascii);

    return 0;
}