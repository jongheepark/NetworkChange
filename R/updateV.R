#' Update layer specific network generation rules
#'
#' Update layer specific network generation rules using optimized matrix operations.
#' Eliminates redundant array permutations for better performance.
#'
#' @param Zb Z - beta.
#' @param U The latent node positions.
#' @param R The dimension of latent space.
#' @param K The dimension of Z.
#' @param s2 The variance of error.
#' @param eV The mean of V.
#' @param iVV The variance of V.
#' @param UTA Indicator of upper triangular array
#'
#' @return A matrix of layer specific network generation rules
#'
#' @export
#'
updateV <- function(Zb, U, R, K, s2, eV, iVV, UTA){
    ## Zero out non-upper triangular elements
    Zb0 <- Zb
    Zb0[!UTA] <- 0

    ## Pre-compute U'U and diagonal correction for Q
    UtU <- crossprod(U)  # t(U) %*% U
    inv_s2 <- 1 / s2

    ## Compute Q = ((U'U)^2 - diag_correction) / 2
    ## This accounts for the symmetric matrix structure
    if(R == 1){
        Q <- (UtU^2 - sum(U^4)) / 2
    } else {
        ## Compute sum of outer products of squared rows
        U_sq <- U^2
        diag_correction <- tcrossprod(colSums(U_sq * U_sq %*% diag(R)))
        ## Simpler: compute directly
        diag_correction <- matrix(0, R, R)
        for(r1 in 1:R){
            for(r2 in 1:R){
                diag_correction[r1, r2] <- sum(U_sq[, r1] * U_sq[, r2])
            }
        }
        Q <- (UtU^2 - diag_correction) / 2
    }

    ## Compute L matrix (K[3] x R) more efficiently
    ## L[t, r] = sum over i < j of Zb[i,j,t] * U[i,r] * U[j,r]
    L <- matrix(0, K[3], R)
    for(t in 1:K[3]){
        Zb_t <- Zb0[,,t]
        for(r in 1:R){
            L[t, r] <- sum(Zb_t * outer(U[,r], U[,r]))
        }
    }

    ## Posterior covariance and mean
    cV <- solve(Q * inv_s2 + iVV)
    ## Prior mean contribution: each row of cE uses same eV
    prior_contrib <- matrix(rep(as.vector(iVV %*% eV), K[3]), nrow = K[3], byrow = TRUE)
    cE <- (L * inv_s2 + prior_contrib) %*% cV
    V <- rmn(cE, diag(K[3]), cV)
    return(V)
}
