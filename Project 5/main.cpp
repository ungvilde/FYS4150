# include "simulate_schrodinger.hpp"

int main() {
    // Read in the input file
    arma::mat system_conditions;
    system_conditions.load("system_conditions.txt", arma::raw_ascii);

    // Define the variables
    int M = system_conditions(0, 0);
    double h = system_conditions(1, 0);
    cx_double r = system_conditions(2, 0);
    cx_double i = system_conditions(3, 0);
    
    // Define the vectors and matrices
    cx_vec a((M - 2)*(M - 2));
    cx_vec b((M - 2)*(M - 2));
    sp_cx_mat A((M - 2)*(M - 2), (M - 2)*(M - 2));
    sp_cx_mat B((M - 2)*(M - 2), (M - 2)*(M - 2));

    // Call the functions
    matrix_setup(r, a, b, M, A, B);
    matrix_setup_v(M, h, dt, V, r, a, b);
    time_loop(A, B, u, b, a);


}
