#' Compare Log Marginal Likelihood
#' 
#' Compare Log Marginal Likelihood
#' 
#' @param outlist List of NetworkChange objects
#' 
#' @return A matrix of log marginal likelihoods. 
#'
#' @seealso \code{\link{WaicCompare}}
#'
#' @references Siddhartha Chib. 1995. ``Marginal Likelihood from the Gibbs Output.''
#' \emph{Journal of the American Statistical Association}. 90: 1313-1321.
#'
#' Jong Hee Park and Yunkyun Sohn. 2019. "Detecting Structural Change
#' in Longitudinal Network Data." \emph{Bayesian Analysis}. Forthcoming.
#'
#' Sumio Watanabe. 2010. "Asymptotic equivalence of Bayes cross validation and widely
#' applicable information criterion in singular learning theory."
#' \emph{Journal of Machine Learning Research}. 11: 3571-3594.


#' @export
MarginalCompare <- function(outlist){
    N.model <- length(outlist)
    marg <- lapply(outlist, attr, "logmarglike")
    breaks <- lapply(outlist, attr, "m")
    outm <- matrix(marg, 1, N.model)
    colnames(outm) <- paste0("break ", breaks) 
    return(outm)
}    
