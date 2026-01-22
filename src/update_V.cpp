// update_V.cpp - V update functions for NetworkChange
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Forward declaration of rmn_cpp from utils.cpp
arma::mat rmn_cpp(const arma::mat& M, const arma::mat& Srow, const arma::mat& Scol);

//' Update Layer-specific Network Generation Rules (C++ version)
//'
//' Updates V using optimized matrix operations.
//'
//' @param Zb 3D array (Z - beta), already zeroed below diagonal
//' @param U Latent node positions matrix (N x R)
//' @param R Dimension of latent space
//' @param K Dimensions of Z (N, N, T)
//' @param s2 Error variance (scalar)
//' @param eV Prior mean vector for V
//' @param iVV Prior precision matrix for V
//' @return Updated V matrix (T x R)
//' @keywords internal
// [[Rcpp::export]]
arma::mat updateV_cpp(const arma::cube& Zb, const arma::mat& U, int R,
                       const Rcpp::IntegerVector& K, double s2,
                       const arma::vec& eV, const arma::mat& iVV) {
    int T_dim = K[2];
    double inv_s2 = 1.0 / s2;

    // Pre-compute U'U
    arma::mat UtU = U.t() * U;

    // Compute diagonal correction for Q
    arma::mat U_sq = arma::square(U);
    arma::mat diag_correction(R, R, arma::fill::zeros);
    for (int r1 = 0; r1 < R; r1++) {
        for (int r2 = 0; r2 < R; r2++) {
            diag_correction(r1, r2) = arma::accu(U_sq.col(r1) % U_sq.col(r2));
        }
    }

    // Q = ((U'U)^2 - diag_correction) / 2
    arma::mat Q = (arma::square(UtU) - diag_correction) / 2.0;

    // Compute L matrix (T_dim x R)
    // L[t, r] = sum over all i,j of Zb[i,j,t] * U[i,r] * U[j,r]
    arma::mat L(T_dim, R, arma::fill::zeros);
    for (int t = 0; t < T_dim; t++) {
        arma::mat Zb_t = Zb.slice(t);
        for (int r = 0; r < R; r++) {
            arma::vec u_r = U.col(r);
            L(t, r) = arma::accu(Zb_t % (u_r * u_r.t()));
        }
    }

    // Posterior covariance
    arma::mat cV = arma::inv_sympd(Q * inv_s2 + iVV);

    // Prior contribution (replicated for each row)
    arma::vec prior_vec = iVV * eV;
    arma::mat prior_contrib(T_dim, R);
    for (int t = 0; t < T_dim; t++) {
        prior_contrib.row(t) = prior_vec.t();
    }

    // Posterior mean
    arma::mat cE = (L * inv_s2 + prior_contrib) * cV;

    // Sample from matrix normal
    return rmn_cpp(cE, arma::eye(T_dim, T_dim), cV);
}

//' Update V from Change-point Network Process (C++ version)
//'
//' Updates regime-specific V matrices using optimized operations.
//'
//' @param ns Number of hidden regimes
//' @param U List of regime-specific latent node positions
//' @param V Current V matrix (not used, kept for interface compatibility)
//' @param Zm List of regime-specific (Z - beta) arrays, zeroed below diagonal
//' @param Km List of regime-specific dimensions
//' @param R Dimension of latent space
//' @param s2 Vector of regime-specific error variances
//' @param eV List of regime-specific prior means for V
//' @param iVV List of regime-specific prior precision matrices for V
//' @return List of regime-specific V matrices
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List updateVm_cpp(int ns, const Rcpp::List& U, const arma::mat& V,
                         const Rcpp::List& Zm, const Rcpp::List& Km, int R,
                         const arma::vec& s2, const Rcpp::List& eV,
                         const Rcpp::List& iVV) {
    Rcpp::List Vm(ns);

    for (int j = 0; j < ns; j++) {
        // Get regime-specific data
        arma::mat Uj = Rcpp::as<arma::mat>(U[j]);
        arma::cube Zj = Rcpp::as<arma::cube>(Zm[j]);
        Rcpp::IntegerVector Km_j = Km[j];
        arma::vec eV_j = Rcpp::as<arma::vec>(eV[j]);
        arma::mat iVV_j = Rcpp::as<arma::mat>(iVV[j]);

        int T_j = Km_j[2];
        double inv_s2_j = 1.0 / s2(j);

        // Pre-compute U'U
        arma::mat UjtUj = Uj.t() * Uj;

        // Compute diagonal correction for Q
        arma::mat Uj_sq = arma::square(Uj);
        arma::mat diag_correction(R, R, arma::fill::zeros);
        for (int r1 = 0; r1 < R; r1++) {
            for (int r2 = 0; r2 < R; r2++) {
                diag_correction(r1, r2) = arma::accu(Uj_sq.col(r1) % Uj_sq.col(r2));
            }
        }

        // Q = ((U'U)^2 - diag_correction) / 2
        arma::mat Q = (arma::square(UjtUj) - diag_correction) / 2.0;

        // Compute L matrix (T_j x R)
        arma::mat L(T_j, R, arma::fill::zeros);
        for (int t = 0; t < T_j; t++) {
            arma::mat Zj_t = Zj.slice(t);
            for (int r = 0; r < R; r++) {
                arma::vec u_r = Uj.col(r);
                L(t, r) = arma::accu(Zj_t % (u_r * u_r.t()));
            }
        }

        // Posterior covariance
        arma::mat cV = arma::inv_sympd(Q * inv_s2_j + iVV_j);

        // Prior contribution
        arma::vec prior_vec = iVV_j * eV_j;
        arma::mat prior_contrib(T_j, R);
        for (int t = 0; t < T_j; t++) {
            prior_contrib.row(t) = prior_vec.t();
        }

        // Posterior mean
        arma::mat cE = (L * inv_s2_j + prior_contrib) * cV;

        // Sample from matrix normal
        Vm[j] = rmn_cpp(cE, arma::eye(T_j, T_j), cV);
    }

    return Vm;
}
