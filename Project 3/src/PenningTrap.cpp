// Definitions for the functions in the PenningTrap class

#include <vector>      // For vector
#include <string>      // For string
#include <stdlib.h>    // For rand (from C). For more powerful random number generation in C++ we should use <random>
#include <stdexcept>   // For runtime_error

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

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r)
{
    // computed as gradient of - V = -V0 / (2*d**2) * (2z**2 - x**2 - y**2)
    arma::vec E = arma::vec(3, arma::fill::zeros);
    r.print("r vector:");
    E(0) = -2.*r(0);
    E(1) = -2.*r(1);
    E(2) = 4.*r(2);

    return -V0/(2*d*d)* E; //electric field from electrodes
}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r)
{
    arma::vec B = arma::vec(3, arma::fill::zeros);
    B(2) = B0;
    return B;
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
    arma::vec F_other_particles = total_force_particles(i);
    return F_external + F_other_particles;
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt)
{

}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt)
{

}

