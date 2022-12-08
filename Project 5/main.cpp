# include "simulate_schrodinger.hpp"

// Compile and run: g++ main.cpp src/simulate_quantum.cpp -I include -std=c++14 -larmadillo -o main Run: main

int main() {
    // Read in the input file
    // arma::mat system_conditions;
    // system_conditions.load("system_conditions.txt", arma::raw_ascii);

    // Define the variables
    int M = 5;  // Number of points along an axis (including the boundaries)
    cx_vec a((M - 2)*(M - 2));  // Diagonal elements of A
    cx_vec b((M - 2)*(M - 2));  // Diagonal elements of B
    sp_cx_mat A = sp_cx_mat((M - 2)*(M - 2), (M - 2)*(M - 2));  // Complex matrix A
    sp_cx_mat B = sp_cx_mat((M - 2)*(M - 2), (M - 2)*(M - 2));  // Complex matrix B
    cx_double r = 1.0;  // Real unit
    cx_double c;  // Imaginary unit

    // double dt; // Time step
    // double h; // Grid spacing (step size)
    // int index(int i, int j, int M); // Wave index function
    // cx_vec u; // Wave function
    // cx_vec a; // Diagonal elements of A
    // cx_vec b; // Off-diagonal elements of A
    // mat V; // Potential matrix
    // vec x; // Grid points  (x-axis)
    // vec y; // Grid points  (y-axis)
    // vec t; // Time points

    //double h = system_conditions(1, 0);
    //cx_double r = system_conditions(2, 0);
    //cx_double i = system_conditions(3, 0);
    
    // Define the vectors and matrices
    //cx_vec a((M - 2)*(M - 2));
    //cx_vec b((M - 2)*(M - 2));
    //sp_cx_mat A((M - 2)*(M - 2), (M - 2)*(M - 2));
    //sp_cx_mat B((M - 2)*(M - 2), (M - 2)*(M - 2));

    // Call the functions
    matrix_setup(M, r, a, b, A, B);
    
    return 0;
}
