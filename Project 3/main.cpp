#include <iostream>    // For cout, endl
#include <string>      // For string

#include "Particle.hpp"
#include "PenningTrap.hpp"
//#include <armadillo>


// 
// Main program
// 
int main()
{
    double q = 1;
    double m = 1;
    arma::vec r1 = arma::vec("20. 0. 20.");
    arma::vec v1 = arma::vec("0. 25. 0.");

    Particle test_particle(q, m, v1, r1);

    arma::vec r2 = arma::vec("25. 25. 0.");
    arma::vec v2 = arma::vec("0. 40. 5.");

    Particle second_test_particle(q, m, v2, r2);

    PenningTrap test_trap(1, 25, 500);

    test_trap.add_particle(test_particle);
    test_trap.add_particle(second_test_particle);
    
    arma::vec F1 = test_trap.total_force_external(0);
    F1.print("External fields force:");

    arma::vec F2 = test_trap.total_force_particles(0);
    F2.print("Other particles force:");

    arma::vec F = test_trap.total_force(0);   
    F.print("Total force:");

    return 0;
}