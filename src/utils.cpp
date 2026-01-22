// utils.cpp - Utility functions for NetworkChange
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Multivariate Normal Random Sampling (C++ version)
//'
//' @param n Number of samples
//' @param mu Mean vector
//' @param Sigma Covariance matrix
//' @return Matrix of n samples (rows) from MVN(mu, Sigma)
//' @keywords internal
// [[Rcpp::export]]
arma::mat rMVNorm_cpp(int n, const arma::vec& mu, const arma::mat& Sigma) {
    int p = mu.n_elem;
    arma::mat E = arma::randn(n, p);
    arma::mat cholSigma = arma::chol(Sigma);  // Upper triangular Cholesky
    arma::mat result = E * cholSigma;
    // Add mean to each row
    for (int i = 0; i < n; i++) {
        result.row(i) += mu.t();
    }
    return result;
}

//' Matrix Normal Random Sampling (C++ version)
//'
//' @param M Mean matrix
//' @param Srow Row covariance matrix
//' @param Scol Column covariance matrix
//' @return Random matrix from Matrix Normal(M, Srow, Scol)
//' @keywords internal
// [[Rcpp::export]]
arma::mat rmn_cpp(const arma::mat& M, const arma::mat& Srow, const arma::mat& Scol) {
    int m = Srow.n_rows;
    int n = Scol.n_rows;

    // Compute square root of Srow via eigen decomposition
    arma::vec eigval_row;
    arma::mat eigvec_row;
    arma::eig_sym(eigval_row, eigvec_row, Srow);
    arma::mat Srow_h = eigvec_row * arma::diagmat(arma::sqrt(arma::abs(eigval_row))) * eigvec_row.t();

    // Compute square root of Scol via eigen decomposition
    arma::vec eigval_col;
    arma::mat eigvec_col;
    arma::eig_sym(eigval_col, eigvec_col, Scol);
    arma::mat Scol_h = eigvec_col * arma::diagmat(arma::sqrt(arma::abs(eigval_col))) * eigvec_col.t();

    // Generate standard normal matrix
    arma::mat Z = arma::randn(m, n);

    // Transform to matrix normal
    return Srow_h * Z * Scol_h + M;
}

//' Gram-Schmidt Orthogonalization (C++ version)
//'
//' @param U Input matrix to orthogonalize
//' @return Orthonormalized matrix
//' @keywords internal
// [[Rcpp::export]]
arma::mat GramSchmidt_cpp(const arma::mat& U) {
    int N = U.n_rows;
    int R = U.n_cols;
    arma::mat Q(N, R);

    // First column
    arma::vec A = U.col(0);
    Q.col(0) = A / std::sqrt(arma::dot(A, A));

    // Subsequent columns
    for (int r = 1; r < R; r++) {
        arma::vec B = U.col(r);
        // Subtract projections onto previous vectors
        for (int k = 0; k < r; k++) {
            B = B - arma::dot(Q.col(k), U.col(r)) * Q.col(k);
        }
        Q.col(r) = B / std::sqrt(arma::dot(B, B));
    }

    return Q;
}

//' Column Normalization (C++ version)
//'
//' Normalize columns of a matrix to have unit Euclidean norm.
//'
//' @param U1 Input matrix
//' @return Matrix with normalized columns
//' @keywords internal
// [[Rcpp::export]]
arma::mat Unormal_cpp(const arma::mat& U1) {
    int R = U1.n_cols;
    arma::mat U2 = U1;

    // Compute column norms
    arma::rowvec su = arma::sqrt(arma::sum(arma::square(U1), 0));

    // Normalize each column
    for (int r = 0; r < R; r++) {
        U2.col(r) = U1.col(r) / su(r);
    }

    return U2;
}

//' Standard Normal Matrix Generator (C++ version)
//'
//' @param m Number of rows
//' @param n Number of columns
//' @return m x n matrix of standard normal random values
//' @keywords internal
// [[Rcpp::export]]
arma::mat rsmn_cpp(int m, int n) {
    return arma::randn(m, n);
}
