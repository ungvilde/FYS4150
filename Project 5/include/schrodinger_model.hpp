# ifndef SIMULATE_SCHRODINGER_HPP
# define SIMULATE_SCHRODINGER_HPP

# include <iostream>
# include <string>
# include <sstream>
# include <cmath>
# include <time.h>
# include <armadillo>

using namespace arma;
using namespace std;


// Integer declarations
int index(int i, int j, int M);

// Function declarations
void matrix_setup(cx_double r, cx_vec a, cx_vec b, int M, sp_cx_mat A, sp_cx_mat B);
void matrix_setup_v(int M, double h, int dt, mat V, cx_double r, cx_vec a, cx_vec b);
void time_loop(sp_cx_mat A, sp_cx_mat B, cx_vec u, cx_vec b, cx_vec a);


# endif
