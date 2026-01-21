#' Update regime-specific variance
#'
#' Update regime-specific variance parameter
#' 
#' @param ns The number of hidden states
#' @param Zm The regime-specific holder of Z - beta
#' @param MU The mean array. 
#' @param c0 Scalar shape parameter
#' @param d0 Scalar scale parameter
#' @param Km Regime-specific dimensions
#'
#' @return A scalar for a regime-specific variance
#'
#' @export
#'
updates2m <- function(ns, Zm, MU, c0, d0, Km){
    s2 <- rep(NA, ns)
    ZEE <- as.list(rep(NA, ns))
    for(j in 1:ns){
        ZEE[[j]] <- Zm[[j]] - MU[[j]]
        EE <- c(ZEE[[j]])    
        s2[j] <- 1/rgamma(1, (c0+prod(Km[[j]]))/2, (d0+ sum(EE^2))/2)
    }  
    return(s2)
}
