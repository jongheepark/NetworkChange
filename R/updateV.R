#' Update layer specific network generation rules
#'
#' Update layer specific network generation rules
#'
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
    Zb0 <- Zb
    Zb0[!UTA] <- 0
    if(R == 1){
        Q <- ((t(U)%*%U)^2 - matrix(sum(U^4), R, R))/2
    }
    if(R > 1){ 
        Q <- ((t(U)%*%U)^2-
                  matrix(apply(apply(U^2,1,function(x){x%*%t(x)}),1,sum),R,R))/2
    }
    UU <- aperm(array( apply(U,1,"*",t(U)) ,dim=c(R,K[1],K[1]) ),c(2,3,1))
    ZbP <- aperm(Zb0,c(3,1,2))
    ZUU <- array(apply(UU,3,function(x){apply(ZbP,1,"*",x)}),
                 dim=c(K[1],K[1],K[3],R))
    L <- apply(ZUU,c(3,4),sum)
    cV <- solve( Q/s2 + iVV)
    cE <- (L/s2 + rep(1,K[3])%*%t(eV)%*%iVV)%*%cV
    V <- rmn(cE, diag(K[3]), cV)
    return(V)
}
