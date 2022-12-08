# include "simulate_schrodinger.hpp"

// Compile and run: g++ main.cpp src/simulate_quantum.cpp -I include -std=c++14 -larmadillo -o main Run: main

int main() {
    // Read in the input file
    // arma::mat system_conditions;
    // system_conditions.load("system_conditions.txt", arma::raw_ascii);

    // Define the variables

    // Problem 2 - Part 1
    int M = 5;  // Number of points along an axis (including the boundaries)
    cx_vec a((M - 2)*(M - 2));  // Diagonal elements of A
    cx_vec b((M - 2)*(M - 2));  // Diagonal elements of B
    sp_cx_mat A = sp_cx_mat((M - 2)*(M - 2), (M - 2)*(M - 2));  // Complex matrix A
    sp_cx_mat B = sp_cx_mat((M - 2)*(M - 2), (M - 2)*(M - 2));  // Complex matrix B
    cx_double r = 1.0;  // Real unit
    cx_double c;  // Imaginary unit

    // Problem 2 - Part 2
    double dt; // Time step
    double h; // Grid spacing (step size)
    mat V; // Potential matrix
    
    // Problem 3
    cx_vec u_0; // Wave function
    cx_vec u_next; // Wave function at the next time step
    int N_t; // Number of time steps

    // Problem 4
    cx_vec x; // Grid points  (x-axis)
    cx_vec y; // Grid points  (y-axis)
    cx_double imag; // Imaginary unit
    // int index(int i, int j, int M); // Wave index function
    // cx_vec u; // Wave function
    // cx_vec a; // Diagonal elements of A
    // cx_vec b; // Off-diagonal elements of A
    

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
    
    // Problem 2 - Part 1
    // matrix_setup(M, r, a, b, A, B);
    // A.raw_print("A:")
    // B.brief_pring("B:")

    // Problem 2 - Part 2 (returns std::out_of_range as of 08.12.22)
    // ab_vector_setup(M, h, dt, V, r, a, b);
    // a.raw_print("a:");
    // b.raw_print("b:");

    // Problem 3
    // solve_unext(A, B, u, b, a, N_t, u_next);
    //u_next.raw_print("u_next:");

    // Problem 4
    // initialize_u(M, x, y, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, u_0, imag);
    
    return 0;
}
