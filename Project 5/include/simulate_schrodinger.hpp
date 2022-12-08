# ifndef SIMULATE_SCHRODINGER_HPP
# define SIMULATE_SCHRODINGER_HPP

# include <iostream>
# include <string.h>
# include <sstream>
# include <cmath>
# include <time.h>
# include <armadillo>

using namespace arma;
using namespace std;

// Global function declarations
int index(int i, int j, int M);
void matrix_setup(cx_double r, cx_vec a, cx_vec b, int M, sp_cx_mat A, sp_cx_mat B);
void output_matrix(sp_cx_mat A, int M, string filename);
// void ab_vector_setup(int M, double h, int dt, mat V, cx_double r, cx_vec a, cx_vec b);
// void solve_unext(sp_cx_mat A, sp_cx_mat B, cx_vec u, cx_vec b, cx_vec a);
// void initialize_u(double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y);


# endif
