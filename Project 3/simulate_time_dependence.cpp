#include <iostream>   
#include <string> 
#include <sstream>
#include <cmath>


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

    // Construct Penning trap
    double B0 = 9.65 * 10;
    double V0 = 2.41*1000000;
    double d = 500;
    PenningTrap penning_trap(B0, V0, d);
    penning_trap.add_particle(p1);

    double f = 0.7;
    double w_V = 2.2;

    penning_trap.set_time_dependent_V(f, w_V);

    double dt = 0.001;
    double time = 0;
    int N_steps = 50000;

    //position data (p1)
    arma::vec x_vals_p1 = arma::vec(N_steps, arma::fill::zeros);
    arma::vec y_vals_p1 = arma::vec(N_steps, arma::fill::zeros);
    arma::vec z_vals_p1 = arma::vec(N_steps, arma::fill::zeros);
    // initial values
    x_vals_p1(0) = p1.position()(0);
    y_vals_p1(0) = p1.position()(1);
    z_vals_p1(0) = p1.position()(2);
    arma::vec r = arma::vec(3, arma::fill::zeros);
    // compute position with Forward Euler
    for(int i = 1; i < N_steps; i++)
    {
        time += dt;
        penning_trap.evolve_RK4(dt, time);
        r = penning_trap.p[0].position();
        x_vals_p1(i) = r(0);
        y_vals_p1(i) = r(1);
        z_vals_p1(i) = r(2);
    }
    std::ostringstream oss;
    oss << "_f_"<< f << "_w_" << w_V << "_RK4" << "_dt_" << dt;
    std::string experiment_info = oss.str();

    x_vals_p1.save("data/x_values_time_dependent_p1" + experiment_info +".txt", arma::raw_ascii);
    y_vals_p1.save("data/y_values_time_dependent_p1" + experiment_info +".txt", arma::raw_ascii);
    z_vals_p1.save("data/z_values_time_dependent_p1" + experiment_info +".txt", arma::raw_ascii);

    return 0;
}