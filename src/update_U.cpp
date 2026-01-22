// update_U.cpp - U update functions for NetworkChange
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Forward declarations from utils.cpp
arma::mat rMVNorm_cpp(int n, const arma::vec& mu, const arma::mat& Sigma);
arma::mat GramSchmidt_cpp(const arma::mat& U);
arma::mat Unormal_cpp(const arma::mat& U1);

//' Update Time-constant Latent Node Positions (C++ version)
//'
//' Updates U using optimized matrix operations with pre-computed invariants.
//'
//' @param K Dimensions of Z (N, N, T)
//' @param U Current latent node positions (N x R)
//' @param V Layer-specific network generation rule (T x R)
//' @param R Dimension of latent space
//' @param Zb 3D array (Z - beta)
//' @param s2 Error variance (scalar)
//' @param eU Prior mean vector for U
//' @param iVU Prior precision matrix for U
//' @return Updated U matrix (N x R)
//' @keywords internal
// [[Rcpp::export]]
arma::mat updateU_cpp(const Rcpp::IntegerVector& K, arma::mat U,
                       const arma::mat& V, int R, const arma::cube& Zb,
                       double s2, const arma::vec& eU, const arma::mat& iVU) {
    int N = K[0];
    int T_dim = K[2];
    double inv_s2 = 1.0 / s2;

    // Pre-compute invariants
    arma::mat VtV = V.t() * V;
    arma::vec iVU_eU = iVU * eU;

    // Random permutation of node indices
    arma::uvec perm = arma::randperm(N);

    for (int idx = 0; idx < N; idx++) {
        int i = perm(idx);

        // Create Ui with row i zeroed
        arma::mat Ui = U;
        Ui.row(i).zeros();

        // Compute U'U for this configuration
        arma::mat UtU = Ui.t() * Ui;

        // Q = (U'U) * (V'V) - Hadamard product
        arma::mat Q = UtU % VtV;

        // Compute L using optimized operations
        // L[r] = sum over j,t of Ui[j,r] * V[t,r] * Zb[i,j,t]
        // First extract zi as N x T matrix
        arma::mat zi(N, T_dim);
        for (int t = 0; t < T_dim; t++) {
            zi.col(t) = Zb.slice(t).row(i).t();
        }

        // L = colSums(Ui * (zi %*% V))
        arma::vec L = arma::sum(Ui % (zi * V), 0).t();

        // Posterior covariance and mean
        arma::mat cV = arma::inv_sympd(Q * inv_s2 + iVU);
        arma::vec cE = cV * (L * inv_s2 + iVU_eU);

        // Sample new U[i,]
        arma::mat sample = rMVNorm_cpp(1, cE, cV);
        U.row(i) = sample.row(0);
    }

    return U;
}

//' Update Regime-specific Latent Node Positions (C++ version)
//'
//' Updates U for each regime using optimized matrix operations.
//'
//' @param ns Number of hidden states
//' @param U List of regime-specific latent node positions
//' @param V Combined V matrix
//' @param R Dimension of latent space
//' @param Zm List of regime-specific (Z - beta) arrays
//' @param Km List of regime-specific dimensions
//' @param ej List of regime indicators (binary vectors)
//' @param s2 Vector of regime-specific error variances
//' @param eU List of regime-specific prior means for U
//' @param iVU List of regime-specific prior precision matrices for U
//' @param UL_Normal Normalization method: "Normal" or "Orthonormal"
//' @return List of updated regime-specific U matrices
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List updateUm_cpp(int ns, Rcpp::List U, const arma::mat& V, int R,
                         const Rcpp::List& Zm, const Rcpp::List& Km,
                         const Rcpp::List& ej, const arma::vec& s2,
                         const Rcpp::List& eU, const Rcpp::List& iVU,
                         const std::string& UL_Normal) {

    for (int j = 0; j < ns; j++) {
        // Get regime-specific data
        arma::mat Uj = Rcpp::as<arma::mat>(U[j]);
        arma::cube Zm_j = Rcpp::as<arma::cube>(Zm[j]);
        Rcpp::IntegerVector Km_j = Km[j];
        Rcpp::IntegerVector ej_j = ej[j];
        arma::vec eU_j = Rcpp::as<arma::vec>(eU[j]);
        arma::mat iVU_j = Rcpp::as<arma::mat>(iVU[j]);

        int N_j = Km_j[0];
        int T_j = Km_j[2];
        double inv_s2_j = 1.0 / s2(j);

        // Extract Vj for this regime (rows where ej == 1)
        std::vector<int> regime_indices;
        for (int t = 0; t < ej_j.size(); t++) {
            if (ej_j[t] == 1) {
                regime_indices.push_back(t);
            }
        }
        arma::mat Vj(T_j, R);
        for (int t = 0; t < T_j; t++) {
            Vj.row(t) = V.row(regime_indices[t]);
        }

        // Pre-compute invariants for this regime
        arma::mat VjtVj = Vj.t() * Vj;
        arma::vec iVU_eU_j = iVU_j * eU_j;

        // Random permutation of node indices
        arma::uvec perm = arma::randperm(N_j);

        for (int idx = 0; idx < N_j; idx++) {
            int i = perm(idx);

            // Create Ui with row i zeroed
            arma::mat Ui = Uj;
            Ui.row(i).zeros();

            // Compute U'U for this configuration
            arma::mat UtU = Ui.t() * Ui;

            // Q = (U'U) * (V'V) - Hadamard product
            arma::mat Q = UtU % VjtVj;

            // Compute L
            // Extract zi as N_j x T_j matrix
            arma::mat zi(N_j, T_j);
            for (int t = 0; t < T_j; t++) {
                zi.col(t) = Zm_j.slice(t).row(i).t();
            }

            // L = colSums(Ui * (zi %*% Vj))
            arma::vec L = arma::sum(Ui % (zi * Vj), 0).t();

            // Posterior covariance and mean
            arma::mat cV = arma::inv_sympd(Q * inv_s2_j + iVU_j);
            arma::vec cE = cV * (L * inv_s2_j + iVU_eU_j);

            // Sample new Uj[i,]
            arma::mat sample = rMVNorm_cpp(1, cE, cV);
            Uj.row(i) = sample.row(0);
        }

        U[j] = Uj;
    }

    // UL normalization
    if (UL_Normal == "Normal") {
        for (int j = 0; j < ns; j++) {
            arma::mat Uj = Rcpp::as<arma::mat>(U[j]);
            U[j] = Unormal_cpp(Uj);
        }
    } else if (UL_Normal == "Orthonormal") {
        for (int j = 0; j < ns; j++) {
            arma::mat Uj = Rcpp::as<arma::mat>(U[j]);
            U[j] = GramSchmidt_cpp(Uj);
        }
    }

    return U;
}
