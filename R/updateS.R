#' Update latent states
#'
#' Update latent states
#'
#' @param iter iteration number
#' @param s the most recent latent states
#' @param V Network generation rules
#' @param m The number of breaks
#' @param Zb Z - b
#' @param Zt  Z stacked by time
#' @param Time The length of time
#' @param MU.state UVU for each state
#' @param P Transition matrix
#' @param s2 error variance
#' @param N.upper.tri The number of upper triangular elements
#' @param random.perturb If \code{random.perturb} = TRUE and a single state observation is found,
#' the latent state is randomly selected by equal weights.
#' 
#' @return A list of vectors containing latent states and their probabilities
#'
#' @export
#'
updateS <- function(iter, s, V, m, Zb, Zt, Time, MU.state, P, s2,
         N.upper.tri, random.perturb){
    ns <- m + 1
    MUt <- list()
    for (t in 1:Time){
        Zt[t, ] <- c(Zb[, , t][upper.tri(Zb[, , t])])
    }
    for(j in 1:ns){
        MUt[[j]] <- matrix(NA,  Time, N.upper.tri)
        for (t in 1:Time){
            MUt[[j]][t, ] <- c((MU.state[[j]][, , t])[upper.tri(MU.state[[j]][, , t])])
        }
    }
    ZMUt <- as.list(rep(NA, ns))## ns by T by upper.tri
    for(j in 1:(m+1)){
        ZMUt[[j]] <- Zt - MUt[[j]]
    }
    ## if(fast){
    ## if(iter == 1) {cat("    Fast state sampling starts!     \n")}
    ## state.out <- ULUstateFastSample(m=m, V, s, s2, P=P, local.type, logistic.tune)
    ## s <- state.out$s
    ## ps <- state.out$ps
    ## } else if(sticky == TRUE){        
    ## if(iter == 1) {cat("    Sticky state sampling starts!     \n")}
    ## state.out <- ULUstateSampleSticky(m=m, ZMUt=ZMUt, s2=s2, P=P)
    ## }else{
    ## if(iter == 1) {cat("    Full state sampling starts!     \n")}
    state.out <- ULUstateSample(m=m, s=s, ZMUt=ZMUt, s2=s2, P=P, random.perturb)
    ## s <- state.out$s
    ## ps <- state.out$ps
    ## }
    return(state.out)
}
