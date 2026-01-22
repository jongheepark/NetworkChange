// tensor_ops.cpp - Tensor operations for NetworkChange
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Create 3D Array from Factor Matrices (C++ version)
//'
//' Computes M(i,j,t) = sum_r U1(i,r) * U2(j,r) * V(t,r)
//' This is the core tensor reconstruction operation used in NetworkChange.
//'
//' @param U1 First factor matrix (N1 x R)
//' @param U2 Second factor matrix (N2 x R), typically same as U1 for undirected networks
//' @param V Third factor matrix (T x R)
//' @return 3D array of dimension (N1 x N2 x T)
//' @keywords internal
// [[Rcpp::export]]
arma::cube M_U_cpp(const arma::mat& U1, const arma::mat& U2, const arma::mat& V) {
    int N1 = U1.n_rows;
    int N2 = U2.n_rows;
    int T_dim = V.n_rows;
    int R = U1.n_cols;

    arma::cube M(N1, N2, T_dim, arma::fill::zeros);

    // Compute M[i,j,t] = sum_r U1[i,r] * U2[j,r] * V[t,r]
    // Optimized by computing outer products for each rank
    for (int r = 0; r < R; r++) {
        // For each time point
        for (int t = 0; t < T_dim; t++) {
            // M(:,:,t) += V[t,r] * U1(:,r) * U2(:,r)'
            M.slice(t) += V(t, r) * U1.col(r) * U2.col(r).t();
        }
    }

    return M;
}

//' Create 3D Array from List of Factor Matrices (C++ version)
//'
//' General version that takes a list of factor matrices.
//' For NetworkChange, typically called with list(U, U, V) for undirected networks.
//'
//' @param U_list List of factor matrices
//' @return Array with dimensions corresponding to the rows of each factor matrix
//' @keywords internal
// [[Rcpp::export]]
arma::cube M_U_list_cpp(const Rcpp::List& U_list) {
    // For the 3D case used in NetworkChange
    // U_list should have 3 elements: U1, U2, V
    if (U_list.size() != 3) {
        Rcpp::stop("M_U_list_cpp currently only supports 3D tensors (list of 3 matrices)");
    }

    arma::mat U1 = Rcpp::as<arma::mat>(U_list[0]);
    arma::mat U2 = Rcpp::as<arma::mat>(U_list[1]);
    arma::mat V = Rcpp::as<arma::mat>(U_list[2]);

    return M_U_cpp(U1, U2, V);
}

//' Compute diagonal correction for V update (C++ version)
//'
//' Computes the diagonal correction term sum_i U_sq(i,r1) * U_sq(i,r2)
//' used in updateV and updateVm functions.
//'
//' @param U Factor matrix
//' @return R x R diagonal correction matrix
//' @keywords internal
// [[Rcpp::export]]
arma::mat compute_diag_correction_cpp(const arma::mat& U) {
    int R = U.n_cols;
    arma::mat U_sq = arma::square(U);
    arma::mat diag_corr(R, R);

    for (int r1 = 0; r1 < R; r1++) {
        for (int r2 = 0; r2 < R; r2++) {
            diag_corr(r1, r2) = arma::accu(U_sq.col(r1) % U_sq.col(r2));
        }
    }

    return diag_corr;
}

//' Compute L matrix for V update (C++ version)
//'
//' Computes L(t, r) = sum over i less than j of Z(i,j,t) * U(i,r) * U(j,r)
//' This is a key bottleneck in updateV.
//'
//' @param Z 3D array (N x N x T), already zeroed below diagonal
//' @param U Factor matrix (N x R)
//' @param T_dim Number of time points
//' @param R Number of latent dimensions
//' @return T_dim x R matrix L
//' @keywords internal
// [[Rcpp::export]]
arma::mat compute_L_matrix_cpp(const arma::cube& Z, const arma::mat& U, int T_dim, int R) {
    arma::mat L(T_dim, R, arma::fill::zeros);

    for (int t = 0; t < T_dim; t++) {
        arma::mat Z_t = Z.slice(t);
        for (int r = 0; r < R; r++) {
            // L[t, r] = sum over all i,j of Z_t[i,j] * U[i,r] * U[j,r]
            // Since Z_t is already zeroed below diagonal, this sums only upper triangle
            arma::vec u_r = U.col(r);
            arma::mat outer_prod = u_r * u_r.t();
            L(t, r) = arma::accu(Z_t % outer_prod);
        }
    }

    return L;
}
