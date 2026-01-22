// hmm_state_sampler.cpp - Hidden Markov Model state sampler for NetworkChange
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Log-sum-exp trick for numerical stability
//'
//' @param x Vector of log values
//' @return log(sum(exp(x)))
//' @keywords internal
inline double log_sum_exp(const arma::vec& x) {
    double max_x = x.max();
    return max_x + std::log(arma::accu(arma::exp(x - max_x)));
}

//' Hidden State Sampler (C++ version)
//'
//' Sample hidden states from hidden Markov multilinear model using
//' forward filtering backward sampling algorithm.
//'
//' @param m Number of breaks
//' @param s Current state vector (not used, kept for interface)
//' @param ZMUt List of (Z - MU) matrices for each state
//' @param s2 Vector of error variances for each state
//' @param P Transition probability matrix
//' @param SOS_random Whether to apply single observation state randomization
//' @return List containing: s (state vector), ps (state probabilities), SOS (flag)
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List ULUstateSample_cpp(int m, const arma::ivec& s, const Rcpp::List& ZMUt,
                               const arma::vec& s2, const arma::mat& P,
                               bool SOS_random) {
    // Get dimensions from first element
    arma::mat ZMUt_0 = Rcpp::as<arma::mat>(ZMUt[0]);
    int T_dim = ZMUt_0.n_rows;
    int N = ZMUt_0.n_cols;  // Number of upper triangular elements
    int ns = m + 1;

    // Compute log densities for all states and time points
    // density_log[j, t] = log probability of observation at time t under state j
    arma::mat density_log(ns, T_dim);
    double log_2pi = std::log(2.0 * M_PI);

    for (int j = 0; j < ns; j++) {
        arma::mat ZMUt_j = Rcpp::as<arma::mat>(ZMUt[j]);
        double log_var = std::log(s2(j));
        double inv_2s2 = 1.0 / (2.0 * s2(j));

        for (int t = 0; t < T_dim; t++) {
            // Sum of squared residuals for time t
            double ss = arma::dot(ZMUt_j.row(t), ZMUt_j.row(t));
            density_log(j, t) = -0.5 * N * (log_2pi + log_var) - ss * inv_2s2;
        }
    }

    // Forward filtering
    arma::mat F(T_dim, ns, arma::fill::zeros);  // Filtered probabilities
    arma::vec pr1(ns, arma::fill::zeros);
    pr1(0) = 1.0;  // Initial probability

    for (int t = 0; t < T_dim; t++) {
        arma::vec pstyt1(ns);
        if (t == 0) {
            pstyt1 = pr1;
        } else {
            pstyt1 = F.row(t-1).t();
            pstyt1 = P.t() * pstyt1;
        }

        // Log-space computation for numerical stability
        arma::vec log_unnorm(ns);
        for (int j = 0; j < ns; j++) {
            log_unnorm(j) = std::log(pstyt1(j) + 1e-300) + density_log(j, t);
        }

        // Normalize in log space then exponentiate
        double log_norm = log_sum_exp(log_unnorm);
        F.row(t) = arma::exp(log_unnorm - log_norm).t();
    }

    // Backward sampling
    arma::ivec s_new(T_dim);
    arma::mat ps(T_dim, ns, arma::fill::zeros);

    // Last time point
    ps.row(T_dim - 1) = F.row(T_dim - 1);
    s_new(T_dim - 1) = ns - 1;  // 0-indexed, so ns-1 is the last state

    for (int t = T_dim - 2; t >= 1; t--) {
        int st = s_new(t + 1);  // Next state (0-indexed)

        // Compute unnormalized backward probability
        arma::vec unnorm_pstyn = F.row(t).t() % P.col(st);
        double sum_unnorm = arma::accu(unnorm_pstyn);

        if (sum_unnorm < 1e-300) {
            // Numerical issue - keep same state
            s_new(t) = s_new(t + 1);
        } else {
            arma::vec pstyn = unnorm_pstyn / sum_unnorm;

            if (st == 0) {
                // Already in first state, must stay
                s_new(t) = 0;
            } else {
                // Sample: either stay in current state or move to previous
                double pone = pstyn(st - 1);
                double u = R::runif(0.0, 1.0);
                s_new(t) = (u < pone) ? (st - 1) : st;
            }
            ps.row(t) = pstyn.t();
        }
    }

    // First time point is always state 0
    s_new(0) = 0;

    // Handle single observation states
    bool SOS = false;
    if (SOS_random) {
        // Count occurrences of each state
        arma::ivec state_counts(ns, arma::fill::zeros);
        for (int t = 0; t < T_dim; t++) {
            state_counts(s_new(t))++;
        }

        // Check if any state has only one observation
        bool has_single_obs = false;
        for (int j = 0; j < ns; j++) {
            if (state_counts(j) == 1) {
                has_single_obs = true;
                break;
            }
        }

        if (has_single_obs) {
            // Random reassignment with equal probability
            arma::vec probs(ns, arma::fill::ones);
            probs /= ns;

            // Generate random states and sort
            arma::ivec random_states(T_dim);
            for (int t = 0; t < T_dim; t++) {
                double u = R::runif(0.0, 1.0);
                double cumsum = 0.0;
                for (int j = 0; j < ns; j++) {
                    cumsum += probs(j);
                    if (u < cumsum) {
                        random_states(t) = j;
                        break;
                    }
                }
            }
            s_new = arma::sort(random_states);

            // Check if all states are represented
            arma::ivec new_counts(ns, arma::fill::zeros);
            for (int t = 0; t < T_dim; t++) {
                new_counts(s_new(t))++;
            }

            bool all_present = true;
            for (int j = 0; j < ns; j++) {
                if (new_counts(j) == 0) {
                    all_present = false;
                    break;
                }
            }

            if (!all_present) {
                // Deterministic assignment
                int per_state = T_dim / ns;
                int remainder = T_dim % ns;
                int idx = 0;
                for (int j = 0; j < ns; j++) {
                    int n_this_state = per_state + (j < remainder ? 1 : 0);
                    for (int k = 0; k < n_this_state && idx < T_dim; k++) {
                        s_new(idx++) = j;
                    }
                }
            }

            SOS = true;
        }
    }

    // Convert to 1-indexed for R
    arma::ivec s_r = s_new + 1;

    return Rcpp::List::create(
        Rcpp::Named("s") = s_r,
        Rcpp::Named("ps") = ps,
        Rcpp::Named("SOS") = SOS
    );
}
