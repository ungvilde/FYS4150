#include "utils.hpp"


// make tridiagonal matrix A = tridag(a, d, a) where A is NxN
arma::mat make_tridiag(int N, double a, double d)
{
  arma::mat A = arma::mat(N, N).fill(0.);
  A.diag() = arma::vec(N).fill(d);
  A.diag(-1) = arma::vec(N-1).fill(a);
  A.diag(1) = arma::vec(N-1).fill(a);
  return A;
}

// Return a string with a double in scientific notation
std::string scientific_format(const double d, const int width, const int prec)
{
  std::stringstream ss;
  ss << std::setw(width) << std::setprecision(prec) << std::scientific << d;
  return ss.str();
}


// Return a string with an armadillo vector in scientific notation
std::string scientific_format(const arma::vec& v, const int width, const int prec)
{
  std::stringstream ss;
  for(double elem : v)
  {
    ss << scientific_format(elem, width, prec);
  }
  return ss.str();
}

//find analytic eigenvalue/vector solution to tridiagonal matrix A = tridiag(a,d,a) where A in NxN
void solve_analytic(int N, double a, double d, arma::vec& eigvals, arma::mat& eigvecs)
{ 
  double pi = 3.14159265358979323846; 

  for(int i = 1; i < N+1; i++)
  {
    // compute eigenvalue i
    eigvals(i-1) = d + 2*a * cos(i*pi / (N+1));

    // compute jth element of eigenvector i
    for(int j=1; j < N+1; j++)
    {
      eigvecs(j-1, i-1) = sin(j*i*pi/(N+1));
    }
  }
  // normalise eigenvectors 
  eigvecs = arma::normalise(eigvecs);
}

