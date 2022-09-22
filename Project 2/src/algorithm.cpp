//
// Definitions for the algorithms declared in clever_algorithms.hpp
//

#include "algorithm.hpp"


// Determine the the max off-diagonal element of a symmetric matrix A
// - Saves the matrix element indicies to k and l 
// - Returns absolute value of A(k,l) as the function return value
double max_offdiag_symmetric(const arma::mat& A, int& k, int &l)
{
  int N = A.n_rows;
  k = N-1;
  l = N-2;
  double max_elem = A(k, l);

  for(int i=0; i < N-1; i++)
  {
    for(int j=i+1; j < N; j++)
    {
      double elem = abs(A(i,j));
      
      if (elem > max_elem)
      {
        k=i;
        l=j;
      }
    }
  }
  return A(k, l);
}

// Performs a single Jacobi rotation, to "rotate away" the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l)
{
  int N = A.n_rows;
  double t;

  // compute tau
  double tau = (A(l,l) - A(k, k))/(2*A(k, l));
  //std::cout << "Tau = " << tau << std::endl; //hmmmm
  
  // compute tan=t, cos=c, sine=s
  if(tau > 0.)
  {
    t = 1 / (tau + sqrt(1+tau*tau));
  }
  else{
    t = -1 / (-tau + sqrt(1+tau*tau));
  }

  double c = 1/sqrt(1+t*t);
  double s = c*t;

  // transform current A matrix
  double A_kk = A(k, k);
  double A_ll = A(l, l); 
  A(k, k) = A_kk*c*c -2*A(k,l)*c*s + A_ll*s*s;
  //std::cout << "A_kk = " << A(k,k) << std::endl;
  A(l, l) = A_ll*c*c + 2*A(k,l)*c*s + A_kk*s*s;
  A(k, l) = 0.;
  A(l, k) = 0.;

  //update for all i not equal to k,l
  for(int i=0; i<N; i++)
  {
    if(i != k && i != l)
    {
      // store new values as intermediary variable
      double A_ik_updated = A(i, k)*c - A(i, l)*s;
      double A_il_updated = A(i, l)*c + A(i, k)*s;

      // then do the update
      A(k, i) = A_ik_updated;
      A(i, k) = A_ik_updated;
      A(i, l) = A_il_updated;
      A(l, i) = A_il_updated;
    }
    
  }
  // Finished updating A

  // Now update overall rotation matrix R (where the eigenvectors are)
  for(int i=0; i<N; i++)
  {
    double R_ik_updated = R(i, k)*c - R(i, l)*s;
    double R_il_updated = R(i,l)*c + R(i,k)*s;

    R(i, k) = R_ik_updated;
    R(i,l) = R_il_updated;
  }

  // Rotation finished
}

// Jacobi method eigensolver:
// - Runs jacobo_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"
//   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
// - Stops if it the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter
void jacobi_eigensolver(arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, const int maxiter, int& iterations, bool& converged)
{
    int N = A.n_rows;
    int k;
    int l;
    arma::mat R = arma::mat(N, N, arma::fill::eye); // rotation matrix

    double max_elem = max_offdiag_symmetric(A, k, l);

    //while(max_elem > eps)
    for(int i=0; i<=maxiter; i++)
    {
        if(max_elem < eps) // if convergence was reached before hitting maxiter
        {
            converged = true;
            iterations = i;
        }

        jacobi_rotate(A, R, k, l);
        max_elem = max_offdiag_symmetric(A, k, l);
        
        if(i==maxiter)
        {
            iterations = maxiter;
        }
    }

    A.print("Eigenvalues");
    R.print("Eigenvectors");
}