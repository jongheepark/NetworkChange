#' Update regime-changing regression parameters
#'
#' Update regime-changing beta
#'
#' @param ns  The number of hidden states
#' @param K  The dimensionality of Z
#' @param s  Latent state vector
#' @param s2 The variance of error
#' @param B0 The prior variance of beta
#' @param p  The rank of X
#' @param ZU Z - ULU
#'
#' @return A vector of regime-changing regression parameters
#'
#' @export
#'
#'
updatebm <- function(ns, K, s, s2, B0, p, ZU){
    cV <- 1/(sum(sapply(1:ns, function(j){
        prod(K[1:2])*sum(s == j)/s2[j]})) +  diag(1/B0, p))
    cE <- cV*sum(unlist(lapply(ZU, sum)))
    bhat <- rnorm(1, cE, sqrt(cV)) 
    return(bhat)
}
