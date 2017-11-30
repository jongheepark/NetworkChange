#' Hidden State Sampler
#'
#' 
#' Sample hidden states from hidden Markov multilinear model
#'
#'
#' 
#' @param m The number of break
#' @param s Latent state vector 
#' @param ZMUt Z - MU
#' @param s2 error variance
#' @param P Transition matrix
#' @param SOS.random single observation state random perturbation
#' 
#' @return A list of a state vector, state probabilities, and SOS.random. 
#'
#' @export
#'

ULUstateSample <- function(m, s, ZMUt, s2, P, SOS.random){
    T <- dim(ZMUt[[1]])[1]
    ns <- m+1
    ## ignore standard error 
    density.log <- as(matrix(unlist(sapply(1:T, function(t){
        lapply(ZMUt, function(x){sum(dnorm(x[t,], 0,
                                           sd = sqrt(s2[s[t]]), log=TRUE))})})), ns, T), "mpfr")
    F   <-  as(matrix(NA, T, m+1), "mpfr")     # storage for the Filtered probabilities
    pr1 <-  as(c(1,rep(0, m)), "mpfr")         # initial probability Pr(s=k|Y0, lambda)
    unnorm.pstyt <- py <- as(rep(NA, m+1), "mpfr")
    for (t in 1:T){
        if(t==1) {
            pstyt1 = pr1
        }else {
            pstyt1 <- F[t-1,]%*%P
        }                
        ## par(mfrow=c(1,2));
        unnorm.pstyt    <- pstyt1*exp(density.log[,t])       
        F[t,]  <-  unnorm.pstyt/sum(unnorm.pstyt) # Pr(st|Yt)
    }
    F <- matrix(as.numeric(F), nrow=T, m+1)
    s      <-  matrix(1, T, 1)   ## holder for state variables
    ps     <-  matrix(NA, T, m+1) ## holder for state probabilities
    ps[T,] <-  F[T,]              ## we know last elements of ps and s
    s[T,1] <-  m+1
    ## t      <-  T-1
    for(t in (T-1):2){
        ## while (t>=1){
        st     <-  s[t+1]
        unnorm.pstyn   <-  F[t,]*P[,st]
        if(sum(unnorm.pstyn) == 0){
            ## if unnorm.pstyn is all zero, what to do?
            cat("F", F[t,]," and P", P[,st]," do not match at t = ", t, "\n")
            s[t]   <- s[t+1]
        } else{
            ## normalize into a prob. density
            pstyn   <-  unnorm.pstyn/sum(unnorm.pstyn)
            if (st==1) {
                s[t]<-1
            }else {
                pone    <-  pstyn[st-1]
                s[t]   <-  ifelse (runif(1) < pone, st-1, st)
            }
            ps[t,] <-  pstyn
            ## probabilities pertaining to a certain state                                                    
        }
        ##     t   <-  t-1                      
    }
    SOS <- FALSE
    if(SOS.random & sum(table(s) == 1)){
        s <- sort(sample(1:ns, T, replace=TRUE, prob=rep(1/ns, ns))) ## 
        ## cat("\n A single observation state is sampled and the latent state sampled randomly.\n")
        if(length(unique(s)) != ns){
            s <- sort(rep(1:ns, length=T))
        }
        SOS <- TRUE
    }
    ## name and report outputs
    out <-  list(s, ps,  SOS)
    names(out)<- c("s","ps", "SOS")
    return(out)
}
