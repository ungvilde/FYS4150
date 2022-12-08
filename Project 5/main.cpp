# include "/home/eventob/non-private-repos/FYS4150/Project 5/include/simulate_schrodinger.hpp"

int main() {
    // Read in the input file
    // arma::mat system_conditions;
    // system_conditions.load("system_conditions.txt", arma::raw_ascii);

    // Global variables

    // int i; // Loop variable
    // int j; // Loop variable
    // int M; // Number of grid points
    // double dt; // Time step
    // double h; // Grid spacing (step size)
    // int index(int i, int j, int M); // Wave index function
    // cx_vec u; // Wave function
    // cx_vec a; // Diagonal elements of A
    // cx_vec b; // Off-diagonal elements of A
    // cx_double r; // Real unit
    // cx_double c; // Imaginary unit
    // sp_cx_mat A; // Matrix A
    // sp_cx_mat B; // Matrix B
    // mat V; // Potential matrix
    // vec x; // Grid points  (x-axis)
    // vec y; // Grid points  (y-axis)
    // vec t; // Time points
    
    // Define the variables
    int M = 5;
    cx_vec a((M - 2)*(M - 2));
    cx_vec b((M - 2)*(M - 2));
    sp_cx_mat A = sp_cx_mat((M - 2)*(M - 2), (M - 2)*(M - 2));
    sp_cx_mat B = sp_cx_mat((M - 2)*(M - 2), (M - 2)*(M - 2));
    cx_double r = 1.0;
    //double h = system_conditions(1, 0);
    //cx_double r = system_conditions(2, 0);
    //cx_double i = system_conditions(3, 0);
    
    // Define the vectors and matrices
    //cx_vec a((M - 2)*(M - 2));
    //cx_vec b((M - 2)*(M - 2));
    //sp_cx_mat A((M - 2)*(M - 2), (M - 2)*(M - 2));
    //sp_cx_mat B((M - 2)*(M - 2), (M - 2)*(M - 2));

    // Call the functions
    matrix_setup(r, a, b, M, A, B);
    output_matrix(A, M, "A.txt");
    
    return 0;
}
