#' Sample a starting value of hidden states
#'
#' Sample a starting value of hidden states
#'
#'
#' @param Z Degree-corrected network array data
#' @param Time The length of time. 
#' @param m The number of breaks
#' @param initial.U Initialized U matrix.
#' @param V Initialized V matrix.
#' @param s2 Initialized error variance
#' @param R The dimensionality of latent space
#'
#' @return A state vector
#'
#' @export

startS <- function(Z, Time, m, initial.U, V, s2, R){
    ns <- m + 1
    K <- dim(Z)
    U.pilot <- matrix(rnorm(K[1]*R), K[1], R)      
    pilot <- UV.lsq(Y = Z, R = R, U=U.pilot, V=V, tol=1e-5)
    ## V <- matrix(apply(attr(pilot,  "Vmat"), 2, mean), K[3], R)
    ## y <<- pilot$V[,which.max(apply(abs(V), 2, sum))]
    mse <- apply(pilot$V, 2, function(x){sum((x - mean(x))^2)})
    rank.order <- order(mse, decreasing = T)
    if(R > 1){
        ## ratio of eigenvalues
        ## 
        Vy <- pilot$V[, rank.order[1]]/pilot$V[, rank.order[2]] ## V[, which.max(mse)]
    } else{
        Vy <- pilot$V
    }
    if(sum(is.na(Vy)) > 0){
        if(which(is.na(Vy))+1 == Time){
            Vy <- ifelse(is.na(Vy), Vy[which(is.na(Vy))-1], Vy)
        } else{
            Vy <- ifelse(is.na(Vy), Vy[which(is.na(Vy))+1], Vy)
        }
    }
    cy <- (Vy - mean(Vy))/sd(Vy)
    mydata <- data.frame(cy = cy)
    b0 <- 0
    B0 <- 1
    change <- MCMCregressChange(cy ~ 1, data=mydata, m = m, b0=b0, B0=B0, sigma.mu=1, sigma.var=1,
                                mcmc=100, burnin=100)
    ## s <- apply(attr(change, "s.store"), 2, median) ## only for MCMCpack 1.4-0
    s <- apply(attr(change, "prob.state"), 1, which.max) ## previous versions for safety
    
    if(length(unique(s)) != m+1){
        s <- sort(rep(1:ns, length=Time)) ## sort(sample(1:ns, size=K[3], replace=TRUE, prob=(rep(1, ns))))
        cat("An equi-distant initial state vector is chosen : ", table(s), "\n")                           
    }         
    return(s)
}
