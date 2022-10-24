// The PenningTrap class

#ifndef __PenningTrap_hpp__  
#define __PenningTrap_hpp__

#include <string>
#include <vector>
#include <armadillo>
#include <cmath>
#include <Particle.hpp> 


class PenningTrap
{

  // Public stuff
  public:

    double B0; // magnetic field strength
    double V0; // applied potential
    double d; // characteristic dimension
    double k_e = 1.38935333 * 100000; // Coloumb constant
    std::vector<Particle> p; // particles in the Penning trap
    bool are_interacting=true; // if there is interaction between particles in Penning trap; true by default
    double freq; // angular frequency of time-dependent electric field
    double amplitude; // amplitude  of time-dependent electric field

    // Constructor
    PenningTrap(double B0_in, double V0_in, double d_in);

    // Add a particle to the trap
    void add_particle(Particle p_in);

    // Add n randomly initialised Ca+ particles to the trap
    void add_particle(int n_particles);
    
    // Determine if particles interact in penning trap or not
    void particle_interaction(bool are_interacting);

    // count number of particles in Penning trap
    int num_particles_in_trap();

    // External electric field at point r=(x,y,z)
    arma::vec external_E_field(arma::vec r);  

    // Time dependent external electric field at r=(x,y,z)
    arma::vec external_E_field(arma::vec r, double time); 

    // set frequency and amplitude of time-dependent voltage
    void set_time_dependent_V(double f, double w_V);
    
    // compute time-dependent voltage
    double time_dependent_V(double time); 
    
    // External magnetic field at point r=(x,y,z)
    arma::vec external_B_field(arma::vec r);  

    // Force on particle_i from particle_j
    arma::vec force_particle(int i, int j);

    // The total force on particle_i from the external fields
    arma::vec total_force_external(int i);

    // The total force on particle_i from the external fields with time-dependent E field
    arma::vec total_force_external(int i, double time);

    // The total force on particle_i from the other particles
    arma::vec total_force_particles(int i);

    // The total force on particle_i from both external fields and other particles
    arma::vec total_force(int i);

    // Total force on particle_i from time-dependent external E field, external B field, and other particles
    arma::vec total_force(int i, double time);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void evolve_RK4(double dt);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order with time-dependent external E field
    void evolve_RK4(double dt, double time);

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_forward_Euler(double dt);

    // Evolve the system one time step (dt) using Forward Euler with time-dependent external E field
    void evolve_forward_Euler(double dt, double time);

};


#endif