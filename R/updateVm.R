#' Update V from a change-point network process
#'
#' Update layer specific network generation rules from a change-point network process.
#' Uses optimized matrix operations and eliminates redundant array permutations.
#'
#' @param ns The number of hidden regimes.
#' @param U The latent node positions.
#' @param V The layer-specific network generation rule.
#' @param Zm The holder of Z - beta.
#' @param Km The dimension of regime-specific Z.
#' @param R The dimension of latent space.
#' @param s2 The variance of error.
#' @param eV The mean of V
#' @param iVV The variance of V
#' @param UTA Indicator of upper triangular array
#'
#' @return A matrix of regime-specific layer specific network generation rules
#'
#' @export
#'
updateVm <- function(ns, U, V, Zm, Km, R, s2, eV, iVV, UTA){
    Vm <- vector("list", ns)

    for(j in 1:ns){
        Uj <- U[[j]]
        Zj <- Zm[[j]]
        Zj[!UTA[[j]]] <- 0

        N_j <- Km[[j]][1]  # Number of nodes
        T_j <- Km[[j]][3]  # Number of time points in this regime
        inv_s2_j <- 1 / s2[j]

        ## Pre-compute U'U
        UjtUj <- crossprod(Uj)

        ## Compute Q = ((U'U)^2 - diag_correction) / 2
        if(R == 1){
            Q <- (UjtUj^2 - sum(Uj^4)) / 2
        } else {
            Uj_sq <- Uj^2
            diag_correction <- matrix(0, R, R)
            for(r1 in 1:R){
                for(r2 in 1:R){
                    diag_correction[r1, r2] <- sum(Uj_sq[, r1] * Uj_sq[, r2])
                }
            }
            Q <- (UjtUj^2 - diag_correction) / 2
        }

        ## Compute L matrix (T_j x R) more efficiently
        ## L[t, r] = sum over i < j of Zj[i,j,t] * Uj[i,r] * Uj[j,r]
        L <- matrix(0, T_j, R)
        for(t in 1:T_j){
            Zj_t <- matrix(Zj[,,t], nrow = N_j, ncol = N_j)  # Ensure matrix form
            for(r in 1:R){
                L[t, r] <- sum(Zj_t * outer(Uj[,r], Uj[,r]))
            }
        }

        ## Posterior covariance and mean
        cV <- solve(Q * inv_s2_j + iVV[[j]])
        prior_contrib <- matrix(rep(as.vector(iVV[[j]] %*% eV[[j]]), T_j), nrow = T_j, byrow = TRUE)
        cE <- (L * inv_s2_j + prior_contrib) %*% cV
        Vm[[j]] <- rmn(cE, diag(T_j), cV)
    }
    return(Vm)
}
