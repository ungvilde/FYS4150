#include "simulate_schrodinger.hpp"

using namespace arma;
using namespace std;

// Problem 2 - Part one: Translate indices (i, j) to 1D index (k)
arma::cx_vec flatten_internal(arma::cx_mat U)
{
    int M = U.n_cols; // M
    arma::cx_vec u = arma::cx_vec((M-2)*(M-2), arma::fill::zeros); // u vector
    int k = 0; // for indexing

    for(int j = 1; j < M-1; j++)
    {
        for(int i = 1; i < M-1; i++)
        {
            u(k) = U(i, j);
            k += 1;
        }
    }
    
    return u;
}

// Problem 2 - Part two: Initialize the structure of matrices A and B
arma::sp_cx_mat get_A_mat(arma::cx_vec a, arma::cx_double r)
{
    int N = a.n_elem; // N = (M-2)*(M-2)
    arma::sp_cx_mat A(N, N);
    int K = sqrt(N); // M-2

    A.diag(-K) -= r;
    A.diag(K) -= r;
    A.diag(-1) -= r;
    A.diag(1) -= r;
    A.diag(0) = a;

    int i;

    for(int k=1; k < K; k++)
    {
        i = k * K - 1; // index for submatrix
        A(i, i+1) = 0; // set zero where boundary cond. are zero
        A(i + 1, i) = 0;
    }
    
    return A;
}

arma::sp_cx_mat get_B_mat(arma::cx_vec b, arma::cx_double r)
{
    int N = b.n_elem; // N = (M-2)*(M-2)
    arma::sp_cx_mat B(N, N);
    int K = sqrt(N); // M-2

    B.diag(0) = b;
    B.diag(-1) += r;
    B.diag(1) += r;
    B.diag(-K) += r;
    B.diag(K) += r;

    int i;

    for(int k=1; k < K; k++)
    {
        i = k * K - 1; // index for submatrix
        B(i, i+1) = 0; // set zero where boundary cond. are zero
        B(i+1, i) = 0;
    }

    return B;
}

// Problem 2 - Part three: Fill the matrices A and B with initial potential V and values of a and b
void matrix_AB_setup(double h, double dt, arma::cx_mat V, arma::sp_cx_mat& A, arma::sp_cx_mat& B)
{
    int M = 1/ h + 1;

    arma::cx_vec a((M-2)*(M-2)); // initialise a and b vectors
    arma::cx_vec b((M-2)*(M-2));

    arma::cx_vec v = flatten_internal(V); // flatten V matrix to get elements wrt 1D k-index, v_ij = v_k

    arma::cx_double r;
    arma::cx_double dt_im;

    r.imag(dt / (2*h*h));
    dt_im.imag(dt); // imaginary number i*dt

    for(int k = 0; k < (M-2)*(M-2); k++)
    {
        a(k) = 1.0 + 4.0 * r + dt_im / 2.0 * v(k);
        b(k) = 1.0 - 4.0 * r - dt_im / 2.0 * v(k);
    }

    A = get_A_mat(a, r);
    B = get_B_mat(b, r);
}

// Problem 4 - Part one: Initialize state u_n as Gaussian wavepacket
arma::cx_double gaussian_wavepacket(double x, double y, double sigma_x, double sigma_y, double p_x, double p_y, double x_c, double y_c)
{
    arma::cx_double p_x_im;
    arma::cx_double p_y_im;

    p_x_im.imag(p_x);
    p_y_im.imag(p_y);

    return exp(
        - (x - x_c)*(x - x_c) / (2.0 * sigma_x*sigma_x)
        - (y - y_c)*(y - y_c) / (2.0 * sigma_y*sigma_y)
        + p_x_im * (x - x_c)
        + p_y_im * (y - y_c)
    );
}

// Problem 4 - Part two: Check boundary conditions and normalize u_n to 1
arma::cx_mat initialise_U(double h, double sigma_x, double sigma_y, double p_x, double p_y, double x_c, double y_c)
{
    int M = 1.0 / h + 1; // num. points, including boundary points
    arma::cx_mat U(M, M, arma::fill::zeros);
    double x;
    double y;

    for(int i=1; i <= M-2; i++)
    {
        for(int j=1; j <= M-2; j++)
        {
            x = h*i; 
            y = h*j;

            U(i, j) = gaussian_wavepacket(x, y, sigma_x, sigma_y, p_x, p_y, x_c, y_c);
        }
    }

    // normalise probabilities
    arma::cx_mat U_conj = arma::conj(U);
    arma::cx_double C = arma::accu(U % U_conj);
    U /= sqrt(C); 

    return U;
}

// Problem 5 - Initialize potential V and potential barrier/wall
arma::cx_mat initialise_V(double h, double v_0, int num_slits)
{
    int M = 1/h + 1; // num points
    int thickness_x = 0.02 / h;
    int slit_sep_y = 0.05 / h;
    int slit_aperture_y = 0.05 / h;
    int center = 0.5 / h;

    arma::cx_mat V(M, M, arma::fill::zeros);

    for(int i=center - thickness_x/2; i <= center + thickness_x/2; i++) // fill wall thickness in x-direction
    {
        for(int j=0; j < M; j++)
        {
            if(num_slits == 2 &&
             (j >= center - slit_sep_y/2 - slit_aperture_y && j <= center - slit_sep_y/2 || // lower slit
             j >= center + slit_sep_y/2 && j <= center + slit_sep_y/2 + slit_aperture_y)) // upper slit
            {
                // std::cout << "2 slit barrier" << std::endl;
                // std::cout << j << std::endl;
            }
            else if(num_slits == 1 && j >= center - slit_aperture_y/2 && j <= center + slit_aperture_y/2) // single slit
            {
                // std::cout << "1 slit barrier" << std::endl;
                // std::cout << j << std::endl;
            }
            else if(num_slits == 3 &&
             (j >= center - slit_aperture_y/2 && j <= center + slit_aperture_y/2 || // mid slit
             j >= center + slit_sep_y + slit_aperture_y/2 && j <= center + slit_sep_y + slit_aperture_y/2 + slit_aperture_y || // upper slit
             j <= center - slit_sep_y - slit_aperture_y/2 && j >= center - slit_sep_y - slit_aperture_y/2 - slit_aperture_y // lower slit
             ))
            {
                // std::cout << "3 slit barrier" << std::endl;
                // std::cout << j << std::endl;
            }
            else
            {
                V(i, j) = v_0; // set potential
            }
        }
    }

    return V;
}

// Problem 3 - Solve the matrix equations for u_n+1
arma::cx_cube simulate(double h, double dt, double T, double x_c, double sigma_x, double p_x, double y_c, double sigma_y, double p_y, double v_0, int num_slits)
{
    
    int N_timesteps = T / dt; // num. timesteps
    int M = 1/h + 1; // num. points
    arma::cx_cube data(M-2, M-2, N_timesteps); // we save each U_n as slice in cube
    arma::cx_mat V = initialise_V(h, v_0, num_slits); 

    V.save("data/potential_3_slit.bin", arma::arma_binary); // maybe save V and plot, to include slit configs in report?

    arma::cx_mat U = initialise_U(h, sigma_x, sigma_y, p_x, p_y, x_c, y_c);
    arma::sp_cx_mat A;
    arma::sp_cx_mat B;
    matrix_AB_setup(h, dt, V, A, B);
    arma::cx_vec u_n = flatten_internal(U);

    data.slice(0) = arma::reshape(u_n, M-2, M-2); // save initial state
    
    for(int n=1; n < N_timesteps; n++) // maybe save initialisation as well
    {
        std::cout << "n = " << n << std::endl;
        arma::cx_vec b = B * u_n;
        u_n = arma::spsolve(A, b);      
        data.slice(n) = arma::reshape(u_n, M-2, M-2);
    }

    return data;
}
