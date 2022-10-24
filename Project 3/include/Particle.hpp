// The Particle class

#ifndef __Particle_hpp__  
#define __Particle_hpp__

#include <string>
#include <armadillo>


class Particle
{
  // All public stuff
  public:
  
    double q; // charge
    double m; // mass
    arma::vec r; // position
    arma::vec v; // velocity

    // Constructor
    Particle(double q_in, double m_in, arma::vec r_in, arma::vec v_in);

    // Method that returns the charge
    double charge();

    // Method that returns the mass
    double mass();

    // Method that returns the position
    arma::vec position();

    // Method that returns the velocity
    arma::vec velocity();

    // Method that updates the velocity
    void update_velocity(arma::vec v_update);

    // Method that updates the position
    void update_position(arma::vec r_update);

};

#endif