#include <iostream>   
#include <string> 
#include <sstream>
#include <cmath>


#include "Particle.hpp"
#include "PenningTrap.hpp"

void two_particle_experiment(PenningTrap penning_trap, Particle p1, Particle p2, bool are_interacting, double dt, double tot_time);
void single_particle_experiment(PenningTrap penning_trap, Particle p, double dt, double tot_time, std::string evolve_method);
arma::vec single_particle_analytic_solution(double time, arma::vec v0, arma::vec r0);

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
    single_particle_experiment(penning_trap, p1, 0.001, 50, "FE");
    single_particle_experiment(penning_trap, p1, 0.001, 50, "RK4");

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

void single_particle_experiment(PenningTrap penning_trap, Particle p, double dt, double tot_time, std::string evolve_method)
{
    penning_trap.add_particle(p);
    int N_steps = tot_time/dt;
    double time;
    arma::vec r_true;
    arma::vec r;
    arma::vec r0 = penning_trap.p[0].position(); // initial position
    arma::vec v0 = penning_trap.p[0].velocity(); // initial velocity

    arma::vec rel_error = arma::vec(N_steps, arma::fill::zeros);

    // compute position with Forward Euler
    for(int i = 1; i < N_steps; i++)
    {
        if(evolve_method == "FE")
        {
            penning_trap.evolve_forward_Euler(dt);
            r = penning_trap.p[0].position();
        }

        if(evolve_method == "RK4")
        {
            penning_trap.evolve_RK4(dt);
            r = penning_trap.p[0].position();
        }
        
        time = i*dt;
        r_true = single_particle_analytic_solution(time, v0, r0);

        rel_error(i) = arma::norm(r - r_true) / arma::norm(r_true);
    }
    std::ostringstream oss;
    oss << "_" << evolve_method << "_dt_" << dt;
    std::string experiment_info = oss.str();
    
    rel_error.save("data/rel_error" + experiment_info + ".txt", arma::raw_ascii);
    // simulate single particle with different dt for 50us
    // - plot relative error of r_i at each time step with all the different dt solutions
    // - do this for both RK4 and FE
    // - estimate error convergence
}

arma::vec single_particle_analytic_solution(double time, arma::vec v0, arma::vec r0)
{
    // should use Particle and PenningTrap objects to get this info?
    double q = 1;
    double m = 40.078; 
    double B0 = 9.65 * 10; 
    double V0 = 2.41*1000000;
    double d = 500; 
    
    arma::vec r = arma::vec(3, arma::fill::zeros);
    
    if(v0(0) != 0 || v0(2) != 0 || r0(1)!=0)
    {
        std::cout << "NB! Analytic solution not valid with the given initial conditions." << std::endl;
        // then break function?
    } 

    double w_z = sqrt(2.*q*V0 / (m * d*d));
    double w0 = q*B0 / m;

    // initial conditions
    double z0 = r0(2);
    double x0 = r0(0);
    double dy0 = v0(1);

    r(2) = z0 * cos(w_z * time); //z

    double w_pluss = (w0 + sqrt(w0*w0 - 2*w_z*w_z)) / 2.;
    double w_minus = (w0 - sqrt(w0*w0 - 2*w_z*w_z)) / 2.;

    double A_pluss = (dy0 + w_minus*x0)/(w_minus - w_pluss);
    double A_minus = -(dy0 + w_pluss*x0)/(w_minus - w_pluss);

    r(0) = A_pluss*cos(-w_pluss*time) + A_minus*cos(-w_minus*time); //x
    r(1) = A_pluss*sin(-w_pluss*time) + A_minus*sin(-w_minus*time); //y

    return r;
}