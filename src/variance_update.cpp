// variance_update.cpp - Variance update functions for NetworkChange
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Update Regime-specific Variance (C++ version)
//'
//' Updates regime-specific variance parameters using inverse gamma sampling.
//'
//' @param ns Number of hidden states
//' @param Zm List of regime-specific (Z - beta) arrays
//' @param MU List of mean arrays
//' @param c0 Shape parameter prior
//' @param d0 Scale parameter prior
//' @param Km List of regime-specific dimensions
//' @return Vector of regime-specific variances
//' @keywords internal
// [[Rcpp::export]]
arma::vec updates2m_cpp(int ns, const Rcpp::List& Zm, const Rcpp::List& MU,
                         double c0, double d0, const Rcpp::List& Km) {
    arma::vec s2(ns);

    for (int j = 0; j < ns; j++) {
        // Get regime-specific data
        arma::cube Zm_j = Rcpp::as<arma::cube>(Zm[j]);
        arma::cube MU_j = Rcpp::as<arma::cube>(MU[j]);
        Rcpp::IntegerVector Km_j = Km[j];

        // Compute residuals: ZEE = Zm - MU
        arma::cube ZEE = Zm_j - MU_j;

        // Sum of squared residuals
        double ss = arma::accu(arma::square(ZEE));

        // Product of dimensions
        double prod_Km = (double)Km_j[0] * Km_j[1] * Km_j[2];

        // Sample from inverse gamma
        // In R: 1/rgamma(1, shape, rate)
        // Shape = (c0 + prod(Km)) / 2
        // Rate = (d0 + sum(EE^2)) / 2
        double shape = (c0 + prod_Km) / 2.0;
        double rate = (d0 + ss) / 2.0;

        // rgamma in R uses shape and rate, but R::rgamma uses shape and scale
        // scale = 1/rate
        s2(j) = 1.0 / R::rgamma(shape, 1.0 / rate);
    }

    return s2;
}
