#' Regime-specific latent node positions
#'
#' Update regime-specific latent node positions using optimized matrix operations.
#' Pre-computes invariant quantities outside the node loop for better performance.
#'
#' @param ns The number of latent states
#' @param U THe latent node positions
#' @param V Layer-specific network generation rule.
#' @param R The dimensionality of latent space
#' @param Zm Regim-specific Z - beta
#' @param Km The dimension of regime-specific Z.
#' @param ej Regime indicator.
#' @param s2 The variance of error.
#' @param eU The regim-specific mean of U.
#' @param iVU The regim-specific variance of U.
#' @param UL.Normal Normalization method for U. "Normal" or "Orthonormal" are supported.
#'
#' @return A matrix of regime-specific latent node positions
#'
#' @export
#'
updateUm <- function(ns, U, V, R, Zm, Km, ej, s2, eU, iVU, UL.Normal){
    for(j in 1:ns){
        Vj <- matrix(V[ej[[j]] == 1, ], nrow = sum(ej[[j]]), ncol = R)

        ## Pre-compute invariants for this regime
        VjtVj <- crossprod(Vj)       # t(Vj) %*% Vj - computed once per regime
        iVU_eU_j <- iVU[[j]] %*% eU[[j]]  # Prior contribution
        inv_s2_j <- 1 / s2[j]        # Avoid repeated division

        N_j <- Km[[j]][1]  # Number of nodes
        T_j <- Km[[j]][3]  # Number of time points in this regime

        for(i in sample(N_j)){
            Ui <- U[[j]]
            Ui[i,] <- 0

            ## Compute U'U for this configuration
            UtU <- crossprod(Ui)

            ## Q = (U'U) * (V'V) - Hadamard product
            Q <- UtU * VjtVj

            ## Compute L using optimized operations
            ## Ensure zi is a matrix even for T_j = 1
            zi <- matrix(Zm[[j]][i,,], nrow = N_j, ncol = T_j)
            L <- colSums(Ui * (zi %*% Vj))  # Vectorized computation

            ## Posterior covariance and mean
            cV <- solve(Q * inv_s2_j + iVU[[j]])
            cE <- cV %*% (L * inv_s2_j + iVU_eU_j)
            U[[j]][i,] <- rMVNorm(1, cE, cV)
        }
    }

    ## UL normalization
    if (UL.Normal == "Normal"){
        for(j in 1:ns){
            U[[j]] <- Unormal(U[[j]])
        }
    } else if(UL.Normal == "Orthonormal"){
        for(j in 1:ns){
            U[[j]] <- GramSchmidt(U[[j]])
        }
    }
    return(U)
}
