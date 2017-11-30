#' Update time-constant regression parameters
#'
#' Update time-constant regression parameters
#'
#' 
#' @param Z Degree corrected response tensor
#' @param MU Mean array
#' @param s2 Error variance
#' @param XtX X'X
#' @param b0 Prior mean of beta
#' @param B0 Prior variance of beta
#' 
#' @return A vector of regression parameters
#'
#' @export
#'
updateb <- function(Z, MU, s2, XtX, b0, B0){
    ZU <- Z - MU
    Xtz <- sum(ZU) 
    cV <- 1/(XtX/s2 +  1/B0)
    cE <- cV*(Xtz/s2 + (1/B0)*b0)  
    bhat <- rnorm(1, cE, sqrt(cV)) 
    return(bhat)
}
