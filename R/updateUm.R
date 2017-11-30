#' Regime-specific latent node positions
#'
#' Update regime-specific latent node positions.
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
        Vj <-  matrix(V[ej[[j]] == 1, ], nrow=sum(ej[[j]]), ncol=R)
        for(i in sample(Km[[j]][1])){
            Ui <- U[[j]]
            Ui[i,] <- 0
            VU <-  aperm(array(apply(Ui,1,"*",t(Vj)), dim=c(R, Km[[j]][3], Km[[j]][1])), c(3,2,1))
            zi <- Zm[[j]][i,,]
            L <-  apply(VU*array(rep(zi,R), dim=c(Km[[j]][1], Km[[j]][3], R)), 3, sum) 
            Q <-  (t(Ui)%*%Ui) * (t(Vj)%*%Vj)
            cV <- solve( Q/s2[j] + iVU[[j]] )
            cE <- cV%*%( L/s2[j] + iVU[[j]]%*%eU[[j]])
            U[[j]][i,] <- rMVNorm(1, cE, cV ) 
        }
    }
    ## UL normalization
    if (UL.Normal == "Normal"){
        for(j in 1:ns){
            U[[j]] <- Unormal(U[[j]])
        }
    }else if(UL.Normal == "Orthonormal"){
        for(j in 1:ns){
            U[[j]] <- GramSchmidt(U[[j]])
        }
    }
    return(U)
}
