#include "/home/eventob/non-private-repos/FYS4150/Project 5/include/simulate_schrodinger.hpp"

using namespace arma;
using namespace std;

// Problem 2 - i,j to k where there is only change in j corresponding to (M - 2) (complete)
int index(int i, int j, int M) {
    int k = i + j*(M - 2);
    return k;
}

// Problem 2 - Part 1: number r and two vectors a, b (length 9) to produce matrices A, B (complete)
// sp_cx_mat is a sparse complex matrix (large matrices with many zeros)
// cx_vec is a complex vector
// cx_double is a complex double
// cx_mat is a complex matrix
// & sign to store values directly
void matrix_setup(int M, cx_double r, cx_vec a, cx_vec b, sp_cx_mat &A, sp_cx_mat &B) {
    // Center diagonal with values from a and b
    for (int i = 0; i <= (M - 3); i++){         // 0, 1, 2
        for (int j = 0; j <= (M - 3); j++){     // 0, 1, 2
            A(index(i, j, M), index(i, j, M)) = a(index(i, j, M));
            B(index(i, j, M), index(i, j, M)) = b(index(i, j, M));
        }
        // Second diagonal r from the center diagonal
        for (int j = 0; j <= (M - 4); j++){     // 0, 1
            A(index(i, j, M), index(i, j + 1, M)) = -r;
            A(index(i, j + 1, M), index(i, j, M)) = -r;
            B(index(i, j, M), index(i, j + 1, M)) = r;
            B(index(i, j + 1, M), index(i, j, M)) = r;
        }
    }
    // First diagonal r from the center diagonal
    for (int i = 0; i <= (M - 4); i++){         // 0, 1
        for (int j = 0; j <= (M - 3); j++){      // 0, 1, 2
            A(index(i + 1, j, M), index(i, j, M)) = -r;
            A(index(i, j, M), index(i + 1, j, M)) = -r;
            B(index(i + 1, j, M), index(i, j, M)) = r;
            B(index(i, j, M), index(i + 1, j, M)) = r;
        }
    }
}

// Problem 2 - Part 2: fill two matrices with inputs M, h, dt, V, r, a, b (incomplete)
void ab_vector_setup(int M, double h, int dt, mat V, cx_double r, cx_vec a, cx_vec b){
    // Indecies is out of bounds, where should we use h? (x_i = i*h, etc.)
    cx_double c = cx_double(dt / 2.);        // Correct implementation of c = (i*dt/2)?
    for (int i = 0; i <= V.n_cols - 1;){     // i = columns
        for (int j = 0; j <= V.n_rows - 1;){ // j = rows
            a(index(i, j, M)) = 1. + 4.*r + c*V(i, j);
            b(index(i, j, M)) = 1. + 4.*r + c*V(i, j);
        }
    }
}

// Problem 3 - Find the next u_n+1 from u_n in a time loop (incomplete, i think)
void solve_unext(sp_cx_mat A, sp_cx_mat B, cx_vec u_0, cx_vec b, cx_vec a, int N_t, cx_vec u_next){
    for (int i = 0; i <= N_t - 1; i++){   // Time loop from 0 to N_t - 1
        cx_vec b = B * u_0;                 // Perform matrix multiplication
        cx_vec u_next = spsolve(A, b);    // Solve the matrix equation for u_{n + 1}
    }
}

// Problem 4 - Set up the initial state u_0ij based on the unnormalised Gaussian wave packet
void initialize_u(int M, cx_vec x, cx_vec y, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, cx_vec u_0, cx_double imag){
    // Raw function, missing input to matrix form
    // Split it into x- and y-components?
    // Normalize the initial state so the sum of all states equals 1
    imag = cx_double(0, 1);
    u_0 = exp(-(pow((x - x_c), 2) / (2 * pow((sigma_x), 2))) - (pow((y - y_c), 2) / (2 * pow(sigma_y, 2))) + (imag * p_x *(x - x_c)) + (imag * p_y *(y - y_c)));
    for (int i= 0; i <= (M - 2); i++){
        for (int j = 0; j <= (M - 2); j++){
            u_0(index(i, j, M)) = u_0(i, j);
        }
    }
}

// Problem 5 - Construct the potential matrix V with barriers for single, double and triple slits
void construct_potential(double v0, int M, mat V){
    thickness = 0.02;
    x_center = 0.5;
    y_center = 0.5;
    slit_aperture = 0.05;
    slit_separation_y = 0.05;
    // How will we differentiate between number of slits?
    for (int i = 0; i < 0.02;){
        // Fill the columns of the matrix with v0
    }
}

// Problem 6 - Once every function works by it self, put them together and store the values in a file
*/
