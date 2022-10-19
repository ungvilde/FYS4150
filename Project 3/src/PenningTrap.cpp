// Definitions for the functions in the PenningTrap class

#include <vector>      // For vector
#include <string>      // For string
#include <stdlib.h>    // For rand (from C). For more powerful random number generation in C++ we should use <random>
#include <stdexcept>   // For runtime_error
#include <cmath>

#include "PenningTrap.hpp"
#include "Particle.hpp"


// Constructor 
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
    B0 = B0_in;
    V0 = V0_in;
    d = d_in;
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in)
{
    PenningTrap::p.push_back(p_in);
}

// Add n randomly initialised particles to the trap
void PenningTrap::add_particle(int n_particles)
{
    for(int i=0; i<n_particles; i++)
    {
        arma::vec r = arma::vec(3).randn() *0.1* d;  
        arma::vec v = arma::vec(3).randn() *0.1* d;
        Particle p(1, 40.078, r, v); // m and q should maybe not be hardcoded?
        add_particle(p);
    }   
}

// Determine if particles interact in penning trap or not
void PenningTrap::particle_interaction(bool are_interacting)
{
    PenningTrap::are_interacting = are_interacting;
}

// count number of particles in Penning trap at a given time
int PenningTrap::num_particles_in_trap()
{
    int N = p.size(); // Num particles added to trap
    arma::vec r;
    int N_trapped = 0;

    for(int i=0; i < N; i++)
    {
        r = p.at(i).position();
        if(arma::norm(r) <= d)
        {
            N_trapped += 1;
        }
    }
    return N_trapped;
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r)
{
    arma::vec E = arma::vec(3, arma::fill::zeros);
    if(arma::norm(r) > d)
    {
        // E field is zero outside trap
    } else
    {
        // computed as gradient of - V = -V0 / (2*d**2) * (2z**2 - x**2 - y**2)
        E(0) = -2.*r(0);
        E(1) = -2.*r(1);
        E(2) = 4.*r(2);
    }

    return -V0/(2*d*d)* E; //electric field from electrodes
}

// External electric field at point r=(x,y,z) with time dependence
arma::vec PenningTrap::external_E_field(arma::vec r, double time)
{
    arma::vec E = arma::vec(3, arma::fill::zeros);
    if(arma::norm(r) > d)
    {
        // E field is zero outside trap
    } else
    {
        // computed as gradient of - V = -V0 / (2*d**2) * (2z**2 - x**2 - y**2)
        E(0) = -2.*r(0);
        E(1) = -2.*r(1);
        E(2) = 4.*r(2);

    }
    
    double Vt = time_dependent_V(time);

    return -Vt/(2*d*d)* E; //electric field from electrodes
}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r)
{
    arma::vec B = arma::vec(3, arma::fill::zeros);
    if(arma::norm(r) > d)
    {
        // B field is zero outside trap
    } else
    {
        B(2) = B0;
    }
    return B;
} 

double PenningTrap::time_dependent_V(double time)
{
    return V0*(1 + amplitude*cos(freq*time));
}

// set frequency and amplitude of time-varying voltage
void PenningTrap::set_time_dependent_V(double f, double w_V)
{
    freq = w_V;
    amplitude = f;
}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j)
{
    double q_i = p.at(i).charge();
    arma::vec r_i = p.at(i).position();
    double q_j = p.at(j).charge();
    arma::vec r_j = p.at(j).position();
    arma::vec E = k_e * q_j * (r_i - r_j) / pow(norm(r_i - r_j), 3);

    // Force on particle i from particle j
    return q_i * E;
}

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i)
{
    arma::vec r_i = p.at(i).position();
    arma::vec v_i = p.at(i).velocity();
    double q_i = p.at(i).charge();
    arma::vec E = external_E_field(r_i);
    arma::vec B = external_B_field(r_i);
    
    // Lorentz force
    arma::vec F = q_i * E + q_i * arma::cross(v_i, B);
    return F;
}

arma::vec PenningTrap::total_force_external(int i, double time)
{
    arma::vec r_i = p.at(i).position();
    arma::vec v_i = p.at(i).velocity();
    double q_i = p.at(i).charge();
    arma::vec E = external_E_field(r_i, time);
    arma::vec B = external_B_field(r_i);
    
    // Lorentz force
    arma::vec F = q_i * E + q_i * arma::cross(v_i, B);
    return F;
}

// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i)
{
    // number of particles in Penning trap
    int N = p.size();
   
    arma::vec F = arma::vec(3, arma::fill::zeros);
    
    //sum over contribution from each particle in trap
    for(int j = 0; j < N; j++)
    {
        if (i == j) 
        {
            continue;
        } 
        F += force_particle(i, j);
    }
    return F;
}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i)
{
    arma::vec F_external = total_force_external(i);
    arma::vec F;
    if(are_interacting)
    {
        arma::vec F_other_particles = total_force_particles(i);
        F = F_external + F_other_particles;
    } else
    {
        F = F_external;   
    }

    return F;
}

// total force on particle_i when external E field is time-dependent
arma::vec PenningTrap::total_force(int i, double time)
{
    arma::vec F;
    arma::vec F_external = total_force_external(i, time);
    if(are_interacting)
    {
        arma::vec F_other_particles = total_force_particles(i);
        F = F_external + F_other_particles;
    } else
    {
        F = F_external;   
    }

    return F;
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt)
{
    int N = p.size();
    std::vector<Particle> p_init = p; // temp. copy of particles in Penning trap

    // compute k_v1 and k_r1
    std::vector<arma::vec> k_v1(N);
    std::vector<arma::vec> k_r1(N);

    for(int i=0; i<N; i++)
    {
        arma::vec v = p.at(i).velocity();
        double m = p.at(i).mass();
        arma::vec F = total_force(i);
        k_v1.at(i) = dt * F / m;
        k_r1.at(i) = dt * v; 
    }

    for(int i=0; i<N; i++)
    {
        arma::vec v = p_init.at(i).velocity();
        arma::vec r = p_init.at(i).position();
        p.at(i).update_position(r + 0.5*k_r1.at(i));
        p.at(i).update_velocity(v + 0.5*k_v1.at(i));
    }

    // compute k_v2 and k_r2
    std::vector<arma::vec> k_v2(N);
    std::vector<arma::vec> k_r2(N);

    for(int i=0; i<N; i++)
    {
        arma::vec v = p.at(i).velocity();
        double m = p.at(i).mass();
        arma::vec F = total_force(i);
        k_v2.at(i) = dt * F / m;
        k_r2.at(i) = dt * v; 
    }

    for(int i=0; i<N; i++)
    {
        arma::vec v = p_init.at(i).velocity();
        arma::vec r = p_init.at(i).position();
        p.at(i).update_position(r + 0.5*k_r2.at(i));
        p.at(i).update_velocity(v + 0.5*k_v2.at(i));
    }

    // compute k_v3 and k_r3
    std::vector<arma::vec> k_v3(N);
    std::vector<arma::vec> k_r3(N);

    for(int i=0; i<N; i++)
    {
        arma::vec v = p.at(i).velocity();
        double m = p.at(i).mass();
        arma::vec F = total_force(i);
        k_v3.at(i) = dt * F / m;
        k_r3.at(i) = dt * v; // Euler-Cromer update
    }

    for(int i=0; i<N; i++)
    {
        arma::vec v = p_init.at(i).velocity();
        arma::vec r = p_init.at(i).position();
        p.at(i).update_position(r + k_r3.at(i));
        p.at(i).update_velocity(v + k_v3.at(i));
    }

    // compute k_v4 and k_r4
    std::vector<arma::vec> k_v4(N);
    std::vector<arma::vec> k_r4(N);

    for(int i=0; i<N; i++)
    {
        arma::vec v = p.at(i).velocity();
        double m = p.at(i).mass();
        arma::vec F = total_force(i);
        k_v4.at(i) = dt * F / m; 
        k_r4.at(i) = dt * v; // Euler-Cromer update
    }

    // do final update of position and velocity
    for(int i=0; i<N; i++)
    {
        arma::vec v = p_init.at(i).velocity();
        arma::vec r = p_init.at(i).position();

        arma::vec v_update = v + 1./6*(k_v1.at(i) + 2.*k_v2.at(i) + 2.*k_v3.at(i) + k_v4.at(i));
        arma::vec r_update = r + 1./6*(k_r1.at(i) + 2.*k_r2.at(i) + 2.*k_r3.at(i) + k_r4.at(i));

        p.at(i).update_position(r_update);
        p.at(i).update_velocity(v_update);
    }
    
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt, double time)
{
    int N = p.size();
    std::vector<Particle> p_init = p; // temp. copy of particles in Penning trap

    // compute k_v1 and k_r1
    std::vector<arma::vec> k_v1(N);
    std::vector<arma::vec> k_r1(N);

    for(int i=0; i<N; i++)
    {
        arma::vec v = p.at(i).velocity();
        double m = p.at(i).mass();
        arma::vec F = total_force(i, time);
        k_v1.at(i) = dt * F / m;
        k_r1.at(i) = dt * v; 
    }

    for(int i=0; i<N; i++)
    {
        arma::vec v = p_init.at(i).velocity();
        arma::vec r = p_init.at(i).position();
        p.at(i).update_position(r + 0.5*k_r1.at(i));
        p.at(i).update_velocity(v + 0.5*k_v1.at(i));
    }

    // compute k_v2 and k_r2
    std::vector<arma::vec> k_v2(N);
    std::vector<arma::vec> k_r2(N);

    for(int i=0; i<N; i++)
    {
        arma::vec v = p.at(i).velocity();
        double m = p.at(i).mass();
        arma::vec F = total_force(i, time + 0.5*dt);
        k_v2.at(i) = dt * F / m;
        k_r2.at(i) = dt * v; 
    }

    for(int i=0; i<N; i++)
    {
        arma::vec v = p_init.at(i).velocity();
        arma::vec r = p_init.at(i).position();
        p.at(i).update_position(r + 0.5*k_r2.at(i));
        p.at(i).update_velocity(v + 0.5*k_v2.at(i));
    }

    // compute k_v3 and k_r3
    std::vector<arma::vec> k_v3(N);
    std::vector<arma::vec> k_r3(N);

    for(int i=0; i<N; i++)
    {
        arma::vec v = p.at(i).velocity();
        double m = p.at(i).mass();
        arma::vec F = total_force(i, time + 0.5*dt);
        k_v3.at(i) = dt * F / m;
        k_r3.at(i) = dt * v; // Euler-Cromer update
    }

    for(int i=0; i<N; i++)
    {
        arma::vec v = p_init.at(i).velocity();
        arma::vec r = p_init.at(i).position();
        p.at(i).update_position(r + k_r3.at(i));
        p.at(i).update_velocity(v + k_v3.at(i));
    }

    // compute k_v4 and k_r4
    std::vector<arma::vec> k_v4(N);
    std::vector<arma::vec> k_r4(N);

    for(int i=0; i<N; i++)
    {
        arma::vec v = p.at(i).velocity();
        double m = p.at(i).mass();
        arma::vec F = total_force(i, time + dt);
        k_v4.at(i) = dt * F / m; 
        k_r4.at(i) = dt * v; // Euler-Cromer update
    }

    // do final update of position and velocity
    for(int i=0; i<N; i++)
    {
        arma::vec v = p_init.at(i).velocity();
        arma::vec r = p_init.at(i).position();

        arma::vec v_update = v + 1./6*(k_v1.at(i) + 2.*k_v2.at(i) + 2.*k_v3.at(i) + k_v4.at(i));
        arma::vec r_update = r + 1./6*(k_r1.at(i) + 2.*k_r2.at(i) + 2.*k_r3.at(i) + k_r4.at(i));

        p.at(i).update_position(r_update);
        p.at(i).update_velocity(v_update);
    }
    
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt)
{
    int N = p.size();
    std::vector<arma::vec> r_updates(N);
    std::vector<arma::vec> v_updates(N);

    for(int i=0; i<N; i++)
    {
        arma::vec v = p.at(i).velocity();
        double m = p.at(i).mass();
        arma::vec F = total_force(i);
        arma::mat v_update = v + dt * F / m;

        arma::vec r = p.at(i).position();
        arma::vec r_update = r + dt * v_update; // Euler-Cromer update

        r_updates.at(i) = r_update;
        v_updates.at(i) = v_update;
    }

    for(int i=0; i<N; i++)
    {
        p.at(i).update_position(r_updates.at(i));
        p.at(i).update_velocity(v_updates.at(i));
    }
    
}

void PenningTrap::evolve_forward_Euler(double dt, double time)
{
    int N = p.size();
    std::vector<arma::vec> r_updates(N);
    std::vector<arma::vec> v_updates(N);

    for(int i=0; i<N; i++)
    {
        arma::vec v = p.at(i).velocity();
        double m = p.at(i).mass();
        arma::vec F = total_force(i, time);
        arma::mat v_update = v + dt * F / m;

        arma::vec r = p.at(i).position();
        arma::vec r_update = r + dt * v_update; // Euler-Cromer update

        r_updates.at(i) = r_update;
        v_updates.at(i) = v_update;
    }

    for(int i=0; i<N; i++)
    {
        p.at(i).update_position(r_updates.at(i));
        p.at(i).update_velocity(v_updates.at(i));
    }
}
