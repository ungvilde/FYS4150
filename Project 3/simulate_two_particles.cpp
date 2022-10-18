#include <iostream>   
#include <string> 
#include <sstream>

#include "Particle.hpp"
#include "PenningTrap.hpp"

void two_particle_experiment(PenningTrap penning_trap, Particle p1, Particle p2, bool are_interacting, double dt, double tot_time);

int main()
{
    // particle 1
    double q = 1;
    double m = 40.078; // atomic mass of Ca+ ion
    arma::vec r1 = arma::vec("20. 0. 20.");
    arma::vec v1 = arma::vec("0. 25. 0.");
    Particle p1(q, m, r1, v1);

    // particle 2
    arma::vec r2 = arma::vec("25. 25. 0.");
    arma::vec v2 = arma::vec("0. 40. 5.");
    Particle p2(q, m, r2, v2);
    
    // Construct Penning trap
    double B0 = 9.65 * 10;
    double V0 = 2.41*1000000;
    double d = 500;
    PenningTrap penning_trap(B0, V0, d);

    two_particle_experiment(penning_trap, p1, p2, true, 0.001, 50);
    two_particle_experiment(penning_trap, p1, p2, false, 0.001, 50);

    // Simulate single particle in trap for 50us

    // plot motion of two particles in trap in xy-plane, with and without interaction
    // - need to implement interaction switch X

    // plot phase space for two particles with and without interaction 

    // make 3D plot of (x, y, z) trajectory with and without interaction

    // values needed:
    // - x, y, z values (p1 and p2)
    // velocity in x and z (p1 and p2)

    // simulate single particle with different dt for 50us
    // - plot relative error of r_i at each time step with all the different dt solutions
    // - do this for both RK4 and FE
    // - estimate error convergence

    // without interaction
    

    return 0;
}

void two_particle_experiment(PenningTrap penning_trap, Particle p1, Particle p2, bool are_interacting, double dt, double tot_time)
{
    // could also include option for method (FE/RK4)
    penning_trap.add_particle(p1);
    penning_trap.add_particle(p2);
    penning_trap.particle_interaction(are_interacting);

    int N_steps = tot_time/dt;

    //position data (p1)
    arma::vec x_vals_p1 = arma::vec(N_steps, arma::fill::zeros);
    arma::vec y_vals_p1 = arma::vec(N_steps, arma::fill::zeros);
    arma::vec z_vals_p1 = arma::vec(N_steps, arma::fill::zeros);
    //velocity data (p1)
    arma::vec v_x_vals_p1 = arma::vec(N_steps, arma::fill::zeros);
    arma::vec v_z_vals_p1 = arma::vec(N_steps, arma::fill::zeros);
    
    //position data (p2)
    arma::vec x_vals_p2 = arma::vec(N_steps, arma::fill::zeros);
    arma::vec y_vals_p2 = arma::vec(N_steps, arma::fill::zeros);
    arma::vec z_vals_p2 = arma::vec(N_steps, arma::fill::zeros);
    //velocity data (p2)
    arma::vec v_x_vals_p2 = arma::vec(N_steps, arma::fill::zeros);
    arma::vec v_z_vals_p2 = arma::vec(N_steps, arma::fill::zeros);

    // initial values
    x_vals_p1(0) = p1.position()(0);
    y_vals_p1(0) = p1.position()(1);
    z_vals_p1(0) = p1.position()(2);

    x_vals_p2(0) = p2.position()(0);
    y_vals_p2(0) = p2.position()(1);
    z_vals_p2(0) = p2.position()(2);

    v_x_vals_p1(0) = p1.velocity()(0);
    v_z_vals_p1(0) = p1.velocity()(2);
    v_x_vals_p2(0) = p2.velocity()(0);
    v_z_vals_p2(0) = p2.velocity()(2);
    
    // compute position with Runge-Kutta 4th order
    for(int i = 1; i < N_steps; i++)
    {
        penning_trap.evolve_RK4(dt);

        // p1 values
        arma::vec r = penning_trap.p[0].position();
        arma::vec v = penning_trap.p[0].velocity();
        x_vals_p1(i) = r(0);
        y_vals_p1(i) = r(1);
        z_vals_p1(i) = r(2);
        v_x_vals_p1(i) = v(0);
        v_z_vals_p1(i) = v(2);

        // p2 values
        r = penning_trap.p[1].position();
        v = penning_trap.p[1].velocity();
        x_vals_p2(i) = r(0);
        y_vals_p2(i) = r(1);
        z_vals_p2(i) = r(2);
        v_x_vals_p2(i) = v(0);
        v_z_vals_p2(i) = v(2);
    }

    std::ostringstream oss;
    oss << "_interaction_" << are_interacting << "_dt_" << dt;
    std::string experiment_info = oss.str();

    x_vals_p1.save("data/x_values_RK4_p1" + experiment_info + ".txt", arma::raw_ascii);
    y_vals_p1.save("data/y_values_RK4_p1" + experiment_info + ".txt", arma::raw_ascii);
    z_vals_p1.save("data/z_values_RK4_p1" + experiment_info + ".txt", arma::raw_ascii);

    x_vals_p2.save("data/x_values_RK4_p2" + experiment_info + ".txt", arma::raw_ascii);
    y_vals_p2.save("data/y_values_RK4_p2" + experiment_info + ".txt", arma::raw_ascii);
    z_vals_p2.save("data/z_values_RK4_p2" + experiment_info + ".txt", arma::raw_ascii);

    v_x_vals_p1.save("data/v_x_values_RK4_p1" + experiment_info + ".txt", arma::raw_ascii);
    v_z_vals_p1.save("data/v_z_values_RK4_p1" + experiment_info + ".txt", arma::raw_ascii);

    v_x_vals_p2.save("data/v_x_values_RK4_p2" + experiment_info + ".txt", arma::raw_ascii);
    v_z_vals_p2.save("data/v_z_values_RK4_p2" + experiment_info + ".txt", arma::raw_ascii);

}
