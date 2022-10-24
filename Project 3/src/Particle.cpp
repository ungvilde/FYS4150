// Definitions for the functions in the Particle class

#include "Particle.hpp"


// Constructor
Particle::Particle(double q_in, double m_in, arma::vec r_in, arma::vec v_in)
{
    q = q_in;
    m = m_in;
    r = r_in;
    v = v_in;
}

// Method that returns the charge
double Particle::charge()
{
    return q;
}

// Method that returns the mass
double Particle::mass()
{
    return m;
}

// Method that returns the position
arma::vec Particle::position()
{
    return r;
}

// Method that returns the velocity
arma::vec Particle::velocity()
{
    return v;
}

// Method that updates the velocity
void Particle::update_velocity(arma::vec v_update)
{
    v = v_update;
}

// Method that updates the velocity
void Particle::update_position(arma::vec r_update)
{
    r = r_update;
}
