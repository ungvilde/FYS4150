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

// // Method that returns a string with info
// std::string Particle::info()
// {
//   std::string info_string = "q = " + std::to_string(q) + ", m = " + std::to_string(m); 
// }