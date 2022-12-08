#include "/home/eventob/non-private-repos/FYS4150/Project 5/include/simulate_schrodinger.hpp"

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
    // for-loop without the initial conditions taken into account
    for (int i = 0; i <= (M - 3); i++){
        for (int j = 0; j < (M - 2); j++){
            A(index(i, j, M), index(i, j, M)) = a(index(i, j, M));
            B(index(i, j, M), index(i, j, M)) = b(index(i, j, M));
            
            /*
            // First diagonal with r
            A(index(i, j, M), index(i, j + 1, M)) = -r;
            A(index(i + 1, j, M), index(i, j, M)) = -r;
            B(index(i, j, M), index(i, j + 1, M)) = r;
            B(index(i + 1, j, M), index(i, j, M)) = r;
            */
            /*
            // Second diagonal with r
            A(index(i, j, M), index(i, j + 3, M)) = -r;
            A(index(i + 3, j, M), index(i, j, M)) = -r;
            B(index(i, j, M), index(i, j + 3, M)) = r;
            B(index(i + 3, j, M), index(i, j, M)) = r;
            */
        }
    }
}

void output_matrix(sp_cx_mat A, int M, string filename) {
    ofstream outfile;
    outfile.open(filename);
    for (int i = 0; i <= (M - 2); i++){
        for (int j = 0; j <= (M - 2); j++){
            outfile << A(index(i, j, M), index(i, j, M)) << " ";
        }
        outfile << endl;
    }
    outfile.close();
}

/*
// Problem 2 - fill two matrices with inputs M, h, dt, V, r, a, b
void ab_vector_setup(int M, double h, int dt, mat V, cx_double r, cx_vec a, cx_vec b){
    cx_double c = cx_double(dt / 2.);   // Correct implementation of value?
    for (int i = 0; i <= (M - 2);){
        for (int j = 0; j <=(M - 2);){
            a = 1. + 4.*r + c*V(i, j);  // Missing a(i, j) specs, what determines length of a?
            b = 1. + 4.*r + c*V(i, j);  // Missing b(i, j) specs, what determines length of b?
        }
    }
}


// Problem 3 - Find the next u_n+1 from u_n in a time loop
void solve_unext(sp_cx_mat A, sp_cx_mat B, cx_vec u, cx_vec b, cx_vec a){
    // Missing a loop
    b = B * u;                   // Perform matrix multiplication
    spsolve(A, u_next, b);       // Solve the matrix equation for u_{n + 1}
    
    return u_next;
}

// Problem 4 - Set up the initial state u_0ij based on the unnormalised Gaussian wave packet
void initialize_u(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y){
    // Raw function, missing input to matrix form
    u_0 = exp(-((x - x_c)**2 / (2 * sigma_x ** 2)) - ((y - y_c)**2 / (2 * sigma_y ** 2)) + (i * p_x *(x - x_c)) + (i * p_y *(y - y_c)));
    
    // Normalize the initial state so the sum of all states equals 1

    return u_0;
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