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

void write_file(std::string filename, std::vector< std::vector<double> > values)
{
    // this is a more general and usable version
    // first we create a file and write the data points to this file
    std::ofstream ofile;
    ofile.open(filename);
    
    // width and precision parameters used to format the output
    int width = 25;
    int precision = 15;
    // Loop over steps
    int M = values[0].size();
    for(int j = 0; j < M; j++)
    {
    ofile << scientific_format(values[0][j], width, precision) //xvalues
        << scientific_format(values[1][j], width, precision) //yvalues
        << std::endl;
    }

    // close the data file
    ofile.close();

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

void compare_solutions(arma::mat eigvecs_analytic, arma::vec eigvals_analytic, arma::mat eigvecs, arma::vec eigvals, double eps)
{
  int N = eigvecs.n_rows;
  arma::vec are_equal(N); // vector containing boolean value when comparing eigenvectors

  bool eigvals_equal = arma::approx_equal(eigvals_analytic, eigvals, "absdiff", eps);
  for(int i=0; i<N; i++)
  {
    // check if eigenvectors are the same (also checks vector when scaled by -1)
    are_equal(i) = arma::approx_equal(eigvecs.col(i), eigvecs_analytic.col(i), "absdiff", eps) || arma::approx_equal(eigvecs_analytic.col(i), -1*eigvecs.col(i), "absdiff", eps);
  }
  bool eigvecs_equal = arma::all(are_equal);

  if(eigvals_equal && eigvecs_equal)
  {
    std::cout << "The eigenvalues and eigenvectors are equal." << std::endl;
  } else{
    std::cout << "The eigenvalues and eigenvectors are NOT equal." << std::endl;
  }
}