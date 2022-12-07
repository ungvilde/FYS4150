#include "simulate_schrodinger.hpp"

using namespace arma;
using namespace std;

// Problem 2 - i,j to k where there is only change in j corresponding to (M - 2)
int index(int i, int j, int M) {
    int k = i + j*(M - 2);
    return k;
}

// Problem 2 - number r and two vectors a, b (length 9) to produce matrices A, B
// sp_cx_mat is a sparse complex matrix (large matrices with many zeros)
// cx_vec is a complex vector
// cx_double is a complex double
// cx_mat is a complex matrix
void matrix_setup(cx_double r, cx_vec a, cx_vec b, int M, sp_cx_mat A, sp_cx_mat B) {
    for (int i = 0; i <= (M - 2);){
        for (int j = 0; j <= (M - 2);){
            A((i, j*(M - 2)), (i, j*(M - 2))) = a((i, j*(M - 2)));
            B(i, j*(M - 2), (i, j*(M - 2))) = b((i, j*(M - 2)));
        }
        for (int j = 0; j <= (M - 2);){
            A((i, j*(M - 2)), (i, j*(M - 2))) = -r;
            B((i, j*(M - 2)), (i, j*(M - 2))) = r;

            A((i, (j - 1)*(M - 2)), (i, j*(M - 2))) = -r;
            B((i, j*(M - 2)), (i, (j - 1)*(M - 2))) = r;
        }

    }

}

// Problem 2 - fill two matrices with inputs M, h, dt, V, r, a, b
void matrix_setup_v(int M, double h, int dt, mat V, cx_double r, cx_vec a, cx_vec b){
    cx_double c = cx_double(dt / 2.);   // Correct implementation of value?
    for (int i = 0; i <= (M - 2);){
        for (int j = 0; j <=(M - 2);){
            a = 1. + 4.*r + c*V(i, j);  // Missing a(i, j) specs, what determines length of a?
            b = 1. + 4.*r + c*V(i, j);  // Missing b(i, j) specs, what determines length of b?
        }
    }
}

// Problem 3 - Find the next u_n+1 from u_n in a time loop
void time_loop(sp_cx_mat A, sp_cx_mat B, cx_vec u, cx_vec b, cx_vec a){
    // Missing a loop
    b = B * u;              // Perform matrix multiplication
    spsolve(A, u, b);       // Solve the matrix equation
}

