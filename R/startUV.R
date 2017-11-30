#' Starting values of U and V
#'
#' Initialize starting values of U and V
#'
#' @param Z Degree-corrected network array data.
#' @param R The dimensionality of latent space.
#' @param K The dimensionality of Z.
#'
#' @return A list of U and V
#'
#' @export
#'
startUV <- function(Z, R, K){
    eig.temp <- eigen(apply(Z,c(1,2),mean))
    d2m <- abs(eig.temp$val)
    U0 <- eig.temp$vec[, order(d2m,decreasing=TRUE) ]
    U <- matrix(U0[, 1:R], nrow=nrow(U0), ncol=R)
    V0 <- rep(1,K[3])%*%t(eig.temp$val)
    V1 <- V0[, order(d2m,decreasing=TRUE) ]
    V <- matrix(V1[, 1:R], K[3], R)
    if(is.complex(V)){
        stop("Principal eigenvalues have a complex number. Computation is halted.")
    }
    out <- list(U, V)
}
