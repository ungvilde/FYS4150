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
arma::cx_vec flatten_internal(arma::cx_mat U);
arma::sp_cx_mat get_A_mat(arma::cx_vec a, arma::cx_double r);
arma::sp_cx_mat get_B_mat(arma::cx_vec b, arma::cx_double r);
void matrix_AB_setup(double h, double dt, arma::cx_mat V, arma::sp_cx_mat& A, arma::sp_cx_mat& B);
arma::cx_double gaussian_wavepacket(double x, double y, double sigma_x, double sigma_y, double p_x, double p_y, double x_c, double y_c);
arma::cx_mat initialise_U(double h, double sigma_x, double sigma_y, double p_x, double p_y, double x_c, double y_c);
arma::cx_mat initialise_V(double h, double v_0, int num_slits);
arma::cx_cube simulate(double h, double dt, double T, double x_c, double sigma_x, double p_x, double y_c, double sigma_y, double p_y, double v_0, int num_slits);


# endif