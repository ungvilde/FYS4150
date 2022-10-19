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

    // particle 2
    arma::vec r2 = arma::vec("500. 500 500.");
    arma::vec v2 = arma::vec("0. 40. 5.");
    Particle p2(q, m, r2, v2);
    
    // Construct Penning trap
    double B0 = 9.65 * 10;
    double V0 = 2.41*1000000;
    double d = 500;
    PenningTrap penning_trap(B0, V0, d);

    penning_trap.add_particle(p1);
    arma::vec B = penning_trap.external_B_field(r2);
    std::cout << "B = \n" << B << std::endl;
    return 0;
}