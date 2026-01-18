#' Hidden State Sampler
#'
#' Sample hidden states from hidden Markov multilinear model.
#' Uses vectorized operations for improved performance.
#'
#' @param m The number of break
#' @param s Latent state vector
#' @param ZMUt Z - MU
#' @param s2 error variance
#' @param P Transition matrix
#' @param SOS.random single observation state random perturbation
#'
#' @return A list of a state vector, state probabilities, and SOS.random.
#'
#' @export
#'

ULUstateSample <- function(m, s, ZMUt, s2, P, SOS.random){
    T <- dim(ZMUt[[1]])[1]
    N <- dim(ZMUt[[1]])[2]  # Number of upper triangular elements
    ns <- m + 1

    ## Vectorized computation of log densities
    ## density.log[j, t] = log probability of observation at time t under state j
    density.log <- matrix(NA, ns, T)
    for(j in 1:ns){
        ## Compute log density for all time points under state j at once
        ## Using vectorized rowSums instead of sapply/lapply
        density.log[j, ] <- -0.5 * N * log(2 * pi * s2[j]) -
                            rowSums(ZMUt[[j]]^2) / (2 * s2[j])
    }

    ## Forward filtering
    F <- matrix(NA, T, ns)   # storage for the Filtered probabilities
    pr1 <- c(1, rep(0, m))   # initial probability Pr(s=k|Y0, lambda)

    for (t in 1:T){
        if(t == 1) {
            pstyt1 <- pr1
        } else {
            pstyt1 <- F[t-1, ] %*% P
        }
        unnorm.pstyt <- pstyt1 * exp(density.log[, t])
        F[t, ] <- unnorm.pstyt / sum(unnorm.pstyt)  # Pr(st|Yt)
    }

    ## Backward sampling
    s <- matrix(1, T, 1)      # holder for state variables
    ps <- matrix(NA, T, ns)   # holder for state probabilities
    ps[T, ] <- F[T, ]         # we know last elements of ps and s
    s[T, 1] <- ns

    for(t in (T-1):2){
        st <- s[t+1]
        unnorm.pstyn <- F[t, ] * P[, st]
        if(sum(unnorm.pstyn) == 0){
            cat("F", F[t, ], " and P", P[, st], " do not match at t = ", t, "\n")
            s[t] <- s[t+1]
        } else {
            pstyn <- unnorm.pstyn / sum(unnorm.pstyn)
            if (st == 1) {
                s[t] <- 1
            } else {
                pone <- pstyn[st-1]
                s[t] <- ifelse(runif(1) < pone, st-1, st)
            }
            ps[t, ] <- pstyn
        }
    }

    ## Handle single observation states
    SOS <- FALSE
    if(SOS.random && any(table(s) == 1)){
        s <- sort(sample(1:ns, T, replace = TRUE, prob = rep(1/ns, ns)))
        if(length(unique(s)) != ns){
            s <- sort(rep(1:ns, length = T))
        }
        SOS <- TRUE
    }

    list(s = s, ps = ps, SOS = SOS)
}

#' Hidden State Sampler with precision
#'
#' Sample hidden states from hidden Markov multilinear model with precision using Rmpfr package.
#' Uses vectorized operations for improved performance.
#'
#' @param m The number of break
#' @param s Latent state vector
#' @param ZMUt Z - MU
#' @param s2 error variance
#' @param P Transition matrix
#' @param SOS.random single observation state random perturbation
#'
#' @return A list of a state vector, state probabilities, and SOS.random.
#'
#' @export
#'

ULUstateSample.mpfr <- function(m, s, ZMUt, s2, P, SOS.random){
    T <- dim(ZMUt[[1]])[1]
    N <- dim(ZMUt[[1]])[2]  # Number of upper triangular elements
    ns <- m + 1

    ## Vectorized computation of log densities
    density.log.num <- matrix(NA, ns, T)
    for(j in 1:ns){
        density.log.num[j, ] <- -0.5 * N * log(2 * pi * s2[j]) -
                                rowSums(ZMUt[[j]]^2) / (2 * s2[j])
    }
    density.log <- as(density.log.num, "mpfr")

    ## Forward filtering with high precision
    F <- as(matrix(NA, T, ns), "mpfr")
    pr1 <- as(c(1, rep(0, m)), "mpfr")

    for (t in 1:T){
        if(t == 1) {
            pstyt1 <- pr1
        } else {
            pstyt1 <- F[t-1, ] %*% P
        }
        unnorm.pstyt <- pstyt1 * exp(density.log[, t])
        F[t, ] <- unnorm.pstyt / sum(unnorm.pstyt)
    }

    ## Convert back to numeric for backward sampling
    F <- matrix(as.numeric(F), nrow = T, ncol = ns)

    ## Backward sampling
    s <- matrix(1, T, 1)
    ps <- matrix(NA, T, ns)
    ps[T, ] <- F[T, ]
    s[T, 1] <- ns

    for(t in (T-1):2){
        st <- s[t+1]
        unnorm.pstyn <- F[t, ] * P[, st]
        if(sum(unnorm.pstyn) == 0){
            cat("F", F[t, ], " and P", P[, st], " do not match at t = ", t, "\n")
            s[t] <- s[t+1]
        } else {
            pstyn <- unnorm.pstyn / sum(unnorm.pstyn)
            if (st == 1) {
                s[t] <- 1
            } else {
                pone <- pstyn[st-1]
                s[t] <- ifelse(runif(1) < pone, st-1, st)
            }
            ps[t, ] <- pstyn
        }
    }

    ## Handle single observation states
    SOS <- FALSE
    if(SOS.random && any(table(s) == 1)){
        s <- sort(sample(1:ns, T, replace = TRUE, prob = rep(1/ns, ns)))
        if(length(unique(s)) != ns){
            s <- sort(rep(1:ns, length = T))
        }
        SOS <- TRUE
    }

    list(s = s, ps = ps, SOS = SOS)
}
