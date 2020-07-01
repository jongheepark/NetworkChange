#' Compare WAIC
#' 
#' Compare Widely Applicable Information Criterion
#' 
#' @param outlist List of NetworkChange objects
#'
#' @return Results of WAIC computation
#' 
#' @seealso \code{\link{MarginalCompare}}
#' 
#' @return A matrix of log marginal likelihoods.
#' 
#' @references 
#' Sumio Watanabe. 2010. "Asymptotic equivalence of Bayes cross validation and widely
#' applicable information criterion in singular learning theory."
#' \emph{Journal of Machine Learning Research}. 11: 3571-3594.
#' 
#' Andrew Gelman, Jessica Hwang, and Aki Vehtari. 2014. "Understanding predictive information
#' criteria for Bayesian models." \emph{Statistics and Computing}. 24(6):997-1016.

#' @export
#' 
WaicCompare <- function(outlist){
    N.model <- length(outlist)
    breaks <- lapply(outlist, attr, "m")
    outl <- lapply(outlist, attr, "Waic.out")
    outm <- matrix(unlist(outl), N.model, 8, byrow=TRUE)
    out <- matrix(outm[,1], 1, N.model)
    colnames(out) <- paste0("break ", breaks) 
    return(out)
}    
