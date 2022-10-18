#include <iostream>    

#include "Particle.hpp"
#include "PenningTrap.hpp"

int main()
{
    // particle 1
    double q = 1;
    double m = 40.078; // atomic mass of Ca+ ion
    arma::vec r1 = arma::vec("20. 0. 20.");
    arma::vec v1 = arma::vec("0. 25. 0.");
    Particle p1(q, m, r1, v1);

    // particle 2
    double q = 1;
    double m = 40.078; 
    arma::vec r2 = arma::vec("25. 25. 0.");
    arma::vec v2 = arma::vec("0. 40. 5.");
    Particle p2(q, m, r1, v1);
    
    // Construct Penning trap
    double B0 = 9.65 * 10;
    double V0 = 2.41*1000000;
    double d = 500;
    PenningTrap penning_trap(B0, V0, d);

    // Simulate single particle in trap for 50us

    // plot motion of two particles in trap in xy-plane, with and without interaction
    // - need to implement interaction switch

    // plot phase space for two particles with and without interaction

    // make 3D plot of (x, y, z) trajectory with and without interaction

    // simulate single particle with different dt for 50us
    // - plot relative error of r_i at each time step with all the different dt solutions
    // - do this for both RK4 and FE
    // - estimate error convergence

    penning_trap.add_particle(p1);
    penning_trap.add_particle(p2);

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
        penning_trap.evolve_forward_Euler(dt);
        arma::vec r = penning_trap.p[0].position();
        x_vals(i) = r(0);
        y_vals(i) = r(1);
        z_vals(i) = r(2);
    }

    x_vals.save("data/x_values_FE.txt", arma::raw_ascii);
    y_vals.save("data/y_values_FE.txt", arma::raw_ascii);
    z_vals.save("data/z_values_FE.txt", arma::raw_ascii);

    PenningTrap test_trap2(B0, V0, d);
    test_trap2.add_particle(p1);

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