# include "simulate_schrodinger.hpp"

int main(){
    
    // Declare variables
    double h = 0.005;
    double dt = 2.5 * pow(10, -5);
    double T = 0.008;

    double x_c = 0.25;
    double sigma_x = 0.05;
    double p_x = 200.0;

    double y_c = 0.5;
    double sigma_y = 0.05;
    double p_y = 0.0;

    double v_0 = 0.0;
    int num_slits = 2;

    arma::cx_cube data;

    // Problem 7 - Without barrier
    data = simulate(h, dt, T, x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0, num_slits);
    data.save("data/problem7partA.bin", arma::arma_binary);

    // Problem 7 - With barrier
    v_0 = pow(10, 10);
    sigma_y = 0.1;
    
    data = simulate(h, dt, T, x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0, num_slits);
    data.save("data/problem7partB.bin", arma::arma_binary);

    // Problem 8 - Changes to the initial state sigma_y with barrier
    T = 0.002;
    sigma_y = 0.2;
    
    data = simulate(h, dt, T, x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0, num_slits);
    data.save("data/problem8.bin", arma::arma_binary);

    // Problem 9 - Probability distribution with one and three slits
    num_slits = 1;
    data = simulate(h, dt, T, x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0, num_slits);
    data.save("data/problem9_1_slit.bin", arma::arma_binary);

    num_slits = 3;
    data = simulate(h, dt, T, x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0, num_slits);
    data.save("data/problem9_3_slits.bin", arma::arma_binary);

    return 0;
}
