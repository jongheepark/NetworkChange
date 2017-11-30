#' Update transition matrix
#'
#' Update transition matrix
#'
#' @param s Latent state vector
#' @param ns The number of hidden states
#' @param P Transition matrix
#' @param A0 Prior of transition matrix
#' 
#' @return A transtion matrix
#'
#' @export
#'
updateP <- function(s, ns, P, A0){
    swit  <-  switchg(s) 
    for (j in 1:ns){
        swit1 <-  A0[j,] + swit[j,]        
        pj <-  rdirichlet.cp(1, swit1)
        P[j,] <-  pj
    }
    return(P)
}
