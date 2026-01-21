#' Update time-constant latent node positions
#'
#' Update time-constant latent node positions using optimized matrix operations.
#' Pre-computes invariant quantities outside the node loop for better performance.
#'
#' @param K The dimensionality of Z
#' @param U The most recent draw of latent node positions
#' @param V Layer-specific network generation rule
#' @param R The dimensionality of latent space
#' @param Zb Z - beta
#' @param s2 error variance
#' @param eU The mean of U
#' @param iVU The variance of U
#'
#' @return A matrix of time-constant latent node positions
#'
#' @export
#'

updateU <- function(K, U, V, R, Zb, s2, eU, iVU){
    ## Pre-compute invariants outside the loop
    VtV <- crossprod(V)  # t(V) %*% V - computed once
    iVU_eU <- iVU %*% eU  # Prior contribution - computed once
    inv_s2 <- 1 / s2      # Avoid repeated division

    for(i in sample(K[1])){
        Ui <- U
        Ui[i,] <- 0

        ## Compute U'U incrementally (subtract row i contribution)
        UtU <- crossprod(Ui)

        ## Q = (U'U) * (V'V) - Hadamard product
        Q <- UtU * VtV

        ## Compute L using optimized operations
        ## L[r] = sum over j,t of Ui[j,r] * V[t,r] * Zb[i,j,t]
        zi <- matrix(Zb[i,,], nrow = K[1], ncol = K[3])  # Ensure matrix form
        L <- colSums(Ui * (zi %*% V))  # Vectorized computation

        ## Posterior covariance and mean
        cV <- solve(Q * inv_s2 + iVU)
        cE <- cV %*% (L * inv_s2 + iVU_eU)
        U[i,] <- rMVNorm(1, cE, cV)
    }
    return(U)
}
