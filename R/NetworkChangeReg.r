#' Changepoint regression analysis of a multilinear tensor model
#'
#' NetworkChangeReg implements Bayesian multiple changepoint regression models to network time series data
#' using a multilinear tensor decomposition method

#' @param Y Reponse tensor 
#' @param X regressor tensor 
#' @param R Dimension of latent space. The default is 2. 
#' @param m Number of change point.
#'  If \code{m = 0} is specified, the result should be the same as \code{NetworkStatic}.
#' @param initial.s The starting value of latent state vector. The default is
#' sampling from equal probabilities for all states. 
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of MCMC iterations after burnin.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' MCMC iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen.  If \code{verbose} is greater than 0 the
#' iteration number, the \eqn{\beta} vector, and the error variance are
#' printed to the screen every \code{verbose}th iteration.
#'
#' 
#' @param UL.Normal Transformation of sampled U. Users can choose "NULL", "Normal" or "Orthonormal."
#' "NULL" is no normalization. "Normal" is the standard normalization.
#' "Orthonormal" is the Gram-Schmidt orthgonalization. Default is "NULL."
#'
#' @param plotUU If \code{plotUU = TRUE} and \code{verbose > 0},
#' then the plot of the latent space will be
#' printed to the screen at every \code{verbose}th iteration.
#' The default is \code{plotUU = FALSE}.
#'
#' @param plotZ If \code{plotZ = TRUE} and \code{verbose > 0},
#' then the plot of the degree-corrected input matrix will be
#' printed to the screen with the sampled mean values at every \code{verbose}th iteration.
#' The default is \code{plotUU = FALSE}.
#'
#' @param constant If \code{constant = TRUE}, constant parameter is sampled
#' and saved in the output as attribute \code{bmat}. Default is \code{constant = FALSE}.
#' 
#' @param b0 The prior mean of \eqn{\beta}. This must be a vector. The default value is 0.
#' @param B0 The prior variance of \eqn{\beta}. This must be a matrix.  The default value is 1.
#' 
#' @param c0 = 0.1
#' @param d0 = 0.1
#' 
#' @param u0 \eqn{u_0/2} is the shape parameter for the inverse
#' Gamma prior on variance parameters for U. The default is 10.

#' @param u1 \eqn{u_1/2} is the scale parameter for the
#' inverse Gamma prior on variance parameters for U.
#' The default is 1.
#'
#' 
#' @param v0 \eqn{v_0/2} is the shape parameter for the inverse
#' Gamma prior on variance parameters for V.
#' The default is 10.
#' 
#' @param v1 \eqn{v_1/2} is the scale parameter for the
#' inverse Gamma prior on variance parameters for V.
#' The default is the time length of Y.
#' 
#' @param a \eqn{a} is the shape1 beta prior for transition probabilities. By default,
#' the expected duration is computed and corresponding a and b values are assigned. The expected
#' duration is the sample period divided by the number of states.
    
#' @param b \eqn{b} is the shape2 beta prior for transition probabilities. By default,
#' the expected duration is computed and corresponding a and b values are assigned. The expected
#' duration is the sample period divided by the number of states.
  
#' @param Waic If \code{Waic = TRUE}, the Watanabe information criterion is computed.
#'
#' @return An mcmc object that contains the posterior sample. This object can
#' be summarized by functions provided by the coda package. The object
#' contains an attribute \code{Waic.out} that contains results of WAIC and the log-marginal
#' likelihood of the model (\code{logmarglike}). The object
#' also contains an attribute \code{prob.state} storage matrix that contains the
#' probability of \eqn{state_i} for each period
#'
#' @seealso \code{\link{NetworkStatic}}
#'
#' @references    Jong Hee Park and Yunkyun Sohn. 2017. "Detecting Structural Change
#' in Network Time Series Data using Bayesian Inference." Working Paper.
#'
#' Peter D. Hoff 2011. "Hierarchical Multilinear Models for Multiway Data."
#' \emph{Computational Statistics \& Data Analysis}. 55: 530-543.
#'
#' Siddhartha Chib. 1998. "Estimation and comparison of multiple change-point models."
#' \emph{Journal of Econometrics}. 86: 221-241.
#'
#' Sumio Watanabe. 2010. "Asymptotic equivalence of Bayes cross validation and widely
#' applicable information criterion in singular learning theory."
#' \emph{Journal of Machine Learning Research}. 11: 3571-3594.

#' Siddhartha Chib. 1995. ``Marginal Likelihood from the Gibbs Output.''
#' \emph{Journal of the American Statistical Association}. 90: 1313-1321.

#' @export
#'
#' @examples
#'
#'    \dontrun{
#'    set.seed(1973)
#'    ## Generate an array (30 by 30 by 40) with block transitions
#'    from 2 blocks to 3 blocks
#'    Y <- MakeBlockNetworkChange(n=10, T=40, type ="split")
#'    G <- 100 ## Small mcmc scans to save time
#' 
#'    ## Fit multiple models for break number detection using Bayesian model comparison
#'    out1 <- NetworkChangeReg(Y, X, R=2, m=1, mcmc=G, burnin=G, verbose=G, Waic=TRUE)
#'
#'    ## plot latent node positions
#'    plotU(out1)
#'  
#'    ## plot layer-specific network generation rules
#'    plotV(out1)
#'    }
X.b.pool <- function(X, bhat)
{  
  p <- dim(X)[length(dim(X)) ]
  tmp <- array(0,dim(X)[-length(dim(X))])
  for(j in 1:p) {tmp <- tmp+X[,,,j]*bhat[j] } 
  tmp
}
###


## time specific Xb calculator
X.b.specific <- function(X, bhat)
{  
  p <- dim(X)[length(dim(X)) ] ## number of X variables
  Time <- dim(X)[length(dim(X)) -1] ## number of time dimensions
  ## dim(X)[-length(dim(X))]: dimension of X excluding p
  tmp.array <- array(0,dim(X)[-length(dim(X))])
  tmp <- matrix(0, dim(X)[1], dim(X)[2])
  for(t in 1:Time){
    for(j in 1:p) {
      tmp <- tmp + X[,,t,j]*bhat[t,j]
    }
    tmp.array[,,t] <- tmp
  }
  tmp.array
}
###
X.b.shrink <- function(X, bhat, state){
    p <- dim(X)[length(dim(X)) ] ## number of X variables
    ## ns <- length(unique(state))
    Time <- dim(X)[length(dim(X)) -1] ## number of time dimensions
    ## end <- c(which(diff(state) == 1), Time)
    ## start <- c(1, which(diff(state) == 1)+1)
    tmp.array <- array(0,dim(X)[-length(dim(X))])
    tmp <- matrix(0, dim(X)[1], dim(X)[2])
    for(t in 1:Time){
        for(j in 1:p) {
            tmp <- tmp + X[,,t,j]*bhat[state[t],j]
        }
        tmp.array[,,t] <- tmp
    }
    tmp.array
}

 mcmc=100; burnin=100; verbose=0; thin  = 1;                         
                             UL.Normal = "Orthonormal";
                             Waic=FALSE; b0 = 0; B0 = 1; c0 = NULL; d0 = NULL;
                             u0 = NULL; u1 = NULL; v0 = NULL; v1 = NULL;
                             a = NULL; b = NULL
        
##################################################################################
NetworkChangeReg <- function(Y, X, R=2, m=1, initial.s = NULL,
                             pooling.mode = c("time.pool", "time.specific", "time.shrink"), 
                             mcmc=100, burnin=100, verbose=0, thin  = 1,                         
                             UL.Normal = "Orthonormal",
                             Waic=FALSE, b0 = 0, B0 = 1, c0 = NULL, d0 = NULL,
                             u0 = NULL, u1 = NULL, v0 = NULL, v1 = NULL,
                             a = NULL, b = NULL){
##################################################################################
    
    ## function call
    ptm <- proc.time()
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)

    ## for future use
    fast = FALSE
    sticky = FALSE
    sequential = FALSE
    local.type = "NULL" ## c("NULL", "linear.trend", "logistic"),
    logistic.tune = 0.5
    random.perturb = TRUE
    
    totiter <- mcmc + burnin
    nstore <- mcmc/thin    
    reduce.mcmc <- nstore
    
    ## changepoint priors and inputs
    ns <- m + 1 # number of states
    ## X <- array(1, dim=c(K, 1))
    p <- dim(X)[4]
    K <- dim(Y)  
    Time <- K[3]

    ## Y to Z transformation
    Z <-  array(qnorm( rank(Y, ties.method="random")/(length(Y)+1) ), dim=K)
    for(k in 1:K[3]) {
        Z[,,k] <-  (Z[,,k] + t(Z[,,k]))/sqrt(2)
    }
    ## tmp <- as.list(rep(NA, K[3]))
    ## nodenames <- dimnames(Y)[[1]]
    
    ## X0 is a upper triangle matrix of X
    X0 <- X ;
    for(k in 1:p) {
        tmp <- X0[,,,k] ; tmp[!UTA] <- 0 ; X0[,,,k] <- tmp
    }
    ## XtX.0 <- apply(X0,c(1,2), function(x){x%*%t(x) } ) 
    XtX.middle <- apply(X0,c(1,2,3), function(x){x%*%t(x) } ) 
    XtX <- matrix(apply(XtX.middle, 1, sum), p, p)
    if(p==1) {
        XtX <- matrix(sum(X0^2), p, p)
    }
    
    ## time specific XtX generator
    XtX.specific <- as.list(rep(NA, Time))
    for(t in 1:Time){
        XtX.middle.specific <- apply(X0[,,t,], c(1,2), function(x){x%*%t(x) } ) ## 16 66 66 array
        XtX.specific[[t]] <- matrix(apply(XtX.middle.specific, 1, sum), p, p)
        if(p==1) {
            XtX.specific[[t]] <- matrix(sum(X0[,,t,]^2), p, p)
        }
    }
    rm(XtX.middle.specific)
    
    ## unique values of Y
    uy <- sort(unique(c(Y)))
        
 
    ## prior for changepoint
    P  <-  NetworkChange:::trans.mat.prior(m=m, n=Time, a = 0.9, b= 0.1)
    A0  <-  NetworkChange:::trans.mat.prior(m=m, n=Time, a = a, b = b)
    nss <- 0
  
    ## if (is.null(initial.V)){
    ## out <- startUV(Z, R, K)
    ## initial.U <- out[[1]]
    ## V <- out[[2]]
    ## MU <- M.U(list(U,U,V))
    
    if(is.null(u0)){
        u0 <- 10
    }
    if(is.null(u1)){
        u1 <- 1 
    }
    ## sigma.mu <- mean(apply(V, 2, mean))
    ## sigma.var <- var(apply(V, 2, mean))
    if(is.null(v0)){
        v0 <- 10
        ## v0 <- 4 + 2 * (sigma.mu^2/sigma.var)
    }
    if(is.null(v1)){
        ## v1 <- 1
        v1 <- K[3]
    }
     
################################
    ## beta set up
    ## initialize beta
################################
    ## pooling
    if(pooling.mode == "time.pool"){
        XM <- NULL;
        for(j in 1:p){
            XM <- cbind(XM, c(X[,,,j]))
        }
        ## vectorized regression
        tmp <- lm( c(Z)~-1+ XM )
        bhat <- tmp$coef
        Zb <-  Z- X.b(X, bhat)
        bhat.mat <- matrix(NA, nstore, p)
         
    }else if(pooling.mode == "time.specific"){
        tmp <- as.list(rep(NA, K[3]))
        bhat <- matrix(NA, Time, p)
        for (t in 1:Time){
            XM <- NULL;
            for(j in 1:p){
                XM <- cbind(XM, c(X[,,t,j]))
            }
            ## vectorized regression
            tmp[[t]] <- lm(c(Z[,,t])~-1+ XM )
            bhat[t,] <- tmp[[t]]$coef
        }
        Zb <- Z - X.b.specific(X, bhat)
        bhat.mat <- matrix(NA, nstore, p*Time)

    }else{
        ## time.shrink
        tmp <- as.list(rep(NA, ns))
        bhat <- matrix(NA, ns, p)
        ## equi-distant state vector for initialization
        state <- sort(sample(1:ns, size=K[3], replace=TRUE, prob=(rep(1, ns))))
        ## median.s <- ceiling(apply(attr(mcmcout, "Smat"), 2, median))
        end <- c(which(diff(state) == 1), Time)
        start <- c(1, which(diff(state) == 1)+1)
        for (j in 1:ns){
            XM <- NULL;
            for(h in 1:p){
                XM <- cbind(XM, c(X[,,start[j]:end[j],h]))
            }
            ## vectorized regression
            tmp[[j]] <- lm(c(Z[,,start[j]:end[j]])~ -1 + XM )
            bhat[j,] <- tmp[[j]]$coef
        } 
        Zb <- Z - X.b.shrink(X, bhat, state)
        bhat.mat <- matrix(NA, nstore, p*ns)
    }
           
    if (is.null(initial.s)){
        s <- state## startS(Z, Time, m, initial.U, V, s2=1, R)
    } else{
        s <- initial.s
    }
    
    ## holder 
    ## Zm is a state-specific holder of ZE = Z - bhat
    ## Zm[[1]] is a subset of Z pertaining to state 1
    ## ZU = Z - ULU
    ## ZY is original Z separated by state
    UTA <- Km <- Zm <- ZY <- ZU <- ej <- U <- MU <- MU.state <- Xm <- Vm <- as.list(rep(NA, ns))
    ps.store <- matrix(0, Time, ns)
    
    ## given the state vector, initialize regime specific U and Vm
    for (j in 1:ns){
        ej[[j]] <- as.numeric(s==j)
        Zm[[j]] <- Zb[,,ej[[j]]==1] 
        tmp <- eigen(apply(Zm[[j]], c(1,2), mean))
        d2m <- abs(tmp$val)
        U0 <- tmp$vec[, order(d2m, decreasing=TRUE) ]
        U[[j]] <- matrix(U0[, 1:R], nrow=nrow(U0), ncol=R)
        Vm[[j]] <- matrix(d2m[1:R], sum(s==j), R, byrow=TRUE)
    }
    V <- Reduce(rbind, Vm)
    
    ## initialize MU and MU.state
    ## MU is regime-specific mean matrix, the length of which depends on regime length
    ## MU.state is a full-length mean matrix for state sampling
    for (j in 1:ns){
        MU[[j]] <-  M.U(list(U[[j]],U[[j]], Vm[[j]]))
        MU.state[[j]] <-  M.U(list(U[[j]],U[[j]],V))
    }
    MUU <- abind(MU)

    ## initialize s2 and d0
    if (is.null(c0)){
        c0 <- 1
    }
    if(is.null(d0)) {
        d0 <- var(as.vector(Z - MU.state[[1]]))
    }
   
    s2 <- 1/rgamma(ns, c0/2, (d0)/2)
    Pmat <- matrix(NA, nstore, ns)
    ## cat("scale prior for sigma2: ", d0, "\n")
    
    ## MCMC holders
    ## outlier <- rep(0, T) ## count the number of times of -Inf
    MU.record <- Umat <- s2mat <- iVU <- eU <- eV <- iVV <- eUmat <- iVUmat <- eVmat <- iVVmat <- as.list(rep(NA, ns))
    for(j in 1:ns){
        s2mat[[j]] <- matrix(NA, nstore)
        Umat[[j]] <- matrix(NA, nstore, K[1]*R)
        eUmat[[j]] <- matrix(NA, nstore, R)
        iVUmat[[j]] <- matrix(NA, nstore, R*R)
        eVmat[[j]]  <- matrix(NA, nstore, R)
        iVVmat[[j]] <- matrix(NA, nstore, R*R)
        MU.record[[j]] <- Y*0
        iVU[[j]] <- diag(R)
        eU[[j]] <- rep(u0, R)
        iVV[[j]] <- diag(R)
        eV[[j]] <- rep(v0, R)
    }
    Vmat <- matrix(NA, nstore, R*K[3])
    Smat <- matrix(NA, nstore, K[3])
    
    ## loglike holder
    N.upper.tri <- K[1]*(K[1]-1)/2
    ## Z.loglike <- matrix(NA, mcmc, K[3])
    ## Z.loglike <- as(matrix(NA, mcmc, K[3]), "mpfr")
    if(Waic){
        Z.loglike.array <- array(NA, dim=c(nstore, N.upper.tri, K[3]))
    }
    logmarglike <- loglike <- logmarglike.upper <- loglike.upper <- NA

    Zt <- matrix(NA,  Time,  N.upper.tri)
    UTAsingle <-  upper.tri(Z[,,1])
    ## UTA array: TRUE for upper triangle
    UTAall <- Z*NA
    for(k in 1:K[3]) {
        UTAall[,,k] <-  upper.tri(Z[,,1] )
    } 
    UTAall <- (UTAall==1)
    
    
    Waic.out <- NA
    SOS <- 0
    
    if(verbose !=0){
        cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n")
        cat("\t NetworkChangeReg Sampler Starts! \n")
        ## cat("\t function called: ")
        ## print(call)
        cat("\t degree normalization: ", degree.normal, "\n")
        cat("\t initial states: ", table(s), "\n")
        cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n")
    }
    
#############################################################
    ## MCMC loop starts!
#############################################################
    for(iter in 1:totiter) {

        ## Step 0. Update Z
        ## update Z
        ## pooling
        if(pooling.mode == "time.pool"){
            EZ <-  X.b(X,bhat) + MUU
        }else if(pooling.mode == "time.specific"){
            EZ <- X.b.specific(X, bhat) + MUU
        }else{
            EZ <- X.b.shrink(X, bhat, state=s) + MUU 
        }
        for(y in sample(uy))
        { 
            lb <- suppressWarnings(max(Z[Y<y & UTAall])) 
            ub <- suppressWarnings(min(Z[Y>y & UTAall]))
            z <- qnorm(runif(sum(Y==y), pnorm( lb-EZ[Y==y] ), pnorm(ub-EZ[Y==y])))
            Z[Y==y] <-  EZ[Y==y] + z
            ## print(y)
        }

        
        
        ## Step 1. update ej, Km, Zm
        ## cat("\n---------------------------------------------- \n ")
        ## cat("Step 1. update ej, Km, Zm \n")
        ## cat("\n---------------------------------------------- \n ")
        for (j in 1:ns){
            ej[[j]] <- as.numeric(s==j)
            Km[[j]] <- dim(Zb[,,ej[[j]]==1])
            
            ## in case state j has only 1 unit, force it to be an array with 1 length
            if(is.na(Km[[j]][3])){
                ZY[[j]] <- array(Z[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1)))
                Zm[[j]] <- array(Zb[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1)))
                Xm[[j]] <- array(X[,,ej[[j]]==1, ], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1), p))
                ## ZU[[j]] <- array(Z[,,ej[[j]]==1] - MU[[j]], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1)))
                 ## Ym[[j]] <- array(Y[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], 1))
            } else{
                ZY[[j]] <- Z[,,ej[[j]]==1]
                Zm[[j]] <- Zb[,,ej[[j]]==1]
                Xm[[j]] <- array(X[,,ej[[j]]==1, ], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1), p))
                ## ZU[[j]] <- Z[,,ej[[j]]==1] - MU[[j]]
                ## Ym[[j]] <- Y[,,ej[[j]]==1]
            }
           
            ## return the right dimension info
            Km[[j]] <- dim(Zm[[j]])

            ## UTA array: TRUE for upper triangle
            UTA[[j]] <- Zm[[j]]*NA
            for(k in 1:Km[[j]][3]) {
                UTA[[j]][,,k] <-  upper.tri(Zm[[j]][,,1])
            } 
            UTA[[j]] <- (UTA[[j]]==1)
        }

        
        ## Step 2. update U
        ## cat("\n---------------------------------------------- \n ")
        ## cat("Step 2. update U \n")
        ## cat("\n---------------------------------------------- \n ")
        U <- NetworkChange::updateUm(ns, U, V, R, Zm, Km, ej, s2, eU, iVU, UL.Normal)
        
        ## Step 3. update V
        ## cat("\n---------------------------------------------- \n ")
        ## cat("Step 3. update V \n")
        ## cat("\n---------------------------------------------- \n ")
        Vm <- NetworkChange::updateVm(ns, U, V, Zm, Km, R, s2, eV, iVV, UTA)
        V <- Reduce(rbind, Vm)

        ## update MU
        for(j in 1:ns){
            ## MU is shorter than MU.state. MU.state is a full length.
            MU[[j]] <- M.U(list(U[[j]],U[[j]],Vm[[j]]))
            MU.state[[j]] <- M.U(list(U[[j]],U[[j]],V))
            ZU[[j]] <- ZY[[j]] - MU[[j]]
        }
        MUU <- abind(MU)
        
        ## Step 4. update s2
        ## cat("\n---------------------------------------------- \n ")
        ## cat("Step 4. update s2 \n")
        ## cat("\n---------------------------------------------- \n ")
        s2 <- NetworkChange::updates2m(ns, Zm, MU, c0, d0, Km)
       
        ## update bhat
        ## pooling
        if(pooling.mode == "time.pool"){
            ZE <- Z - MUU
            Xtz <- t(apply(X0,4,c))%*%c(ZE)
            cV <- solve( XtX + diag(1/100,p))
            cE <- cV%*%Xtz
            bhat <- rmvnorm(1,cE,cV)
            Zb <-  Z- X.b(X, bhat)            
        }else if(pooling.mode == "time.specific"){
            for (t in 1:Time){
                ZE <- Z[,,t] - MUU[,,t]
                Xtz <-  t(apply(X0[,,t,],3,c))%*%c(ZE)
                cV <- solve(XtX.specific[[t]] + diag(1/B0, p))
                cE <- cV%*%Xtz
                bhat[t, ] <- rmvnorm(1,cE,cV)
            }
            Zb <- Z - X.b.specific(X, bhat)            
        }else{
            end <- c(which(diff(s) == 1), Time)
            start <- c(1, which(diff(s) == 1)+1)
            for (j in 1:ns){
                ZE <- ZY[[j]] - MU[[j]]
                Xtz <-  t(apply(Xm[[j]],4,c))%*%c(ZE)
                XtX.middle <- apply(Xm[[j]],c(1,2,3), function(x){x%*%t(x) } ) ## 16 66 66  8 array
                XtX <- matrix(apply(XtX.middle, 1, sum), p, p)
                cV <- solve( XtX + diag(1/B0,p))
                cE <- cV%*%Xtz
                ## vectorized regression
                bhat[j,] <- rmvnorm(1,cE,cV)
            } 
            Zb <- Z - X.b.shrink(X, bhat, state=s)
        }
        
        
        ## update hierarchical parameters
        ## hierarchical parameters for U
        for(j in 1:ns){
            SS <-  t(U[[j]]) %*% U[[j]]## (Km[[j]][1]-1)*cov(U[[j]]) + Km[[j]][1]*msi/(Km[[j]][1]+1)
            for(r in 1:R){
                iVU[[j]][r,r] <- 1/rgamma(1, (u0 + K[1])/2, (u1+ SS[r,r])/2)
            }
            eU[[j]] <- c(NetworkChange:::rMVNorm(1,apply(U[[j]],2,sum)/(Km[[j]][1]+1), solve(iVU[[j]])/(Km[[j]][1]+1)))
        }
        
        ## hierarchical parameters for V
        ## V for state j only
        for(j in 1:ns){
            Vs <- matrix(Vm[[j]], nrow=sum(ej[[j]]), ncol=R)
            SS <-  t(Vs)%*%Vs
            for(r in 1:R){
                iVV[[j]][r,r] <- 1/rgamma(1, (v0 + Km[[j]][3])/2, (v1 + SS[r,r])/2)
            }
            eV[[j]] <- c(NetworkChange:::rMVNorm(1,apply(Vs, 2, sum)/(Km[[j]][3]+1),
                                                 solve(iVV[[j]])/(Km[[j]][3]+1)))      
        }
        
        
        ## Step 5. update s
        ## cat("\n---------------------------------------------- \n ")
        ## cat("Step 5. update s \n")
        ## cat("\n---------------------------------------------- \n ")
        state.out <- NetworkChange::updateS(iter, s, V, m, Zb, Zt, Time, MU.state, P, s2,
                             N.upper.tri, random.perturb)
        ## state.out <- updateS(iter, s, V, m, Zb, Zt, Time, fast,
        ##                       MU.state, P, s2, local.type, logistic.tune, N.upper.tri, sticky)
        s <- state.out$s
        ps <- state.out$ps
 
        ## double check 
        if(length(table(s)) < ns){
            ## print(table(s))
            ## cat("Sampled s does not have all states. \n")
            s <- sort(sample(1:ns, size=K[3], replace=TRUE, prob=(rep(1, ns))))
        }
    
        
        ## Step 6. update P
        ## cat("\n---------------------------------------------- \n ")
        ## cat("Step 6. update P \n")
        ## cat("\n---------------------------------------------- \n ")       
        P <- NetworkChange::updateP(s, ns, P, A0)

        
        ## report
        if (verbose!= 0 &iter %% verbose == 0){
            cat("\n----------------------------------------------",'\n')
            cat("    iteration = ", iter, '\n')
            ## cat("    SOS = ", SOS, '\n')
            cat("    beta = ", bhat,'\n')
            if(plotZ == TRUE & plotUU == TRUE){
                if(ns < 4){
                    par(mfrow=c(1, ns+1))
                } else{
                    par(mfrow=c(2, ceiling((ns+1)/2)))
                }
            }
            if(plotZ == TRUE){
                plot(density(c(Z)), lwd=2, main="Density of Z and MU.state")
                for(j in 1:ns){lines(density(c(MU.state[[j]])), col=1+j)}
                legend("topright", paste0("Regime", 1:ns), col=2:(ns+1), lty=1, lwd=1)
            }
            
            if(plotZ == FALSE & plotUU == TRUE){
                par(mfrow=c(1, ns))
            }
            
            for(j in 1:ns){
                cat("    state ", j, "has :", sum(s==j),'\n')
                cat("    sigma2 at state", j, "=", s2[j] ,'\n')
                if(plotUU == TRUE){
                    plot(U[[j]][,1], U[[j]][, 2], pch=19, cex=1); abline(v=0, col=2); abline(h=0, col=2)
                }           
            }
            cat("----------------------------------------------",'\n')
        }
        
        ## save
        if (iter > burnin & (iter-burnin)%%thin == 0){
            nss <- nss + 1

            bhat.mat[iter-burnin, ] <- bhat

            for(j in 1:ns){
                MU.record[[j]] <- MU.record[[j]] + MU.state[[j]]
                s2mat[[j]][(iter-burnin)/thin] <- s2[j]
                Umat[[j]][(iter-burnin)/thin, ] <- as.vector(U[[j]])
                eUmat[[j]][(iter-burnin)/thin, ] <- as.vector(eU[[j]])
                iVUmat[[j]][(iter-burnin)/thin, ] <- as.vector(iVU[[j]])
                eVmat[[j]][(iter-burnin)/thin, ] <- as.vector(eV[[j]])
                iVVmat[[j]][(iter-burnin)/thin, ] <- as.vector(iVV[[j]])
            }
            Vmat[(iter-burnin)/thin, ] <- as.vector(V)
            Smat[(iter-burnin)/thin, ] <- s
            Pmat[(iter-burnin)/thin, ] <- diag(P)
            ps.store <- ps.store + ps
            if(Waic){
                d <- sapply(1:K[3], function(t){dnorm(c(Zb[,,t][UTAsingle]),
                                                      mean = c(MU.state[[s[t]]][,,t][UTAsingle]),
                                                      sd=sqrt(s2[[s[t]]]), log=TRUE)})
                Z.loglike.array[(iter-burnin)/thin, ,] <- d
            }
        }
        
    }## end of MCMC loop
  


    
##########################################################################
    ## Model Diagnostics
##########################################################################
    s.st <- ceiling(apply(Smat, 2, median))
    M.st <-  lapply(MU.record, function(x){x/nss})## MU.record[[j]]/nss      
    if(time.pool){
        bhat.st <- apply(bhat.mat, 2, mean) ## mean((trueY - M.st)^2 )
        Zb.st <- Z - bhat.st 
    }else if(time.specific){
        bhat.st <- apply(bhat.mat, 2, mean)
        Zb.st <- Z - bhat.st 
    }else{
        bhat.st <- apply(bhat.mat, 2, mean)
        Zb.st <- Z - bhat.st 
    }
    
    P.st <- apply(Pmat, 2, mean)
    Km.st <- Km
    for(j in 1:ns){
        Km.st[[j]][3] <- sum(s.st == j)
    }
     
    ## regime specific stars
    U.st <- eU.st <- eV.st <- iVU.st <- iVV.st <- Vm.st <- 
        MU.st <- MU.state.st <- EE.st <- ZU.st <- as.list(rep(NA, ns))
    s2.st <- rep(NA, ns)
    for(j in 1:ns){
        U.st[[j]] <- matrix(apply(Umat[[j]], 2, mean), K[1], R)
        Vm.st[[j]] <- matrix(apply(matrix(Vmat[, s.st == j], nstore, R*sum(s.st == j)), 2, mean), sum(s.st == j), R)
        eU.st[[j]] <-apply(eUmat[[j]], 2, mean)
        eV.st[[j]] <- apply(eVmat[[j]], 2, mean)
        iVU.st[[j]] <- matrix(apply(iVUmat[[j]], 2, mean), R, R)
        iVV.st[[j]] <- matrix(apply(iVVmat[[j]], 2, mean), R, R)
        
        MU.st[[j]] <- M.U(list(U.st[[j]], U.st[[j]], Vm.st[[j]]))
        MU.state.st[[j]] <- M.U(list(U.st[[j]],U.st[[j]], matrix(apply(Vmat, 2, mean), K[3], R)))
        EE.st[[j]] <- c(array(Zb.st[,,s.st == j], dim=Km.st[[j]]) - MU.st[[j]]) ## M.U(list(U,U,V)) 
    }
    V.st <- Reduce(rbind, Vm.st)
    ## apply(V.st, 2, sum)
    ## t(U.st[[1]][,1]) %*% U.st[[1]][,2]
    
    for(j in 1:ns){
        ej[[j]]  <-  as.numeric(s.st==j)
        N.regime <- length(which(s.st == j))    
        s2.st[j] <- mean(s2mat[[j]]) ## mean((trueY - M.st)^2 )
        
        ## in case state j has only 1 unit, force it to be an array with 1 length
        if(is.na(Km.st[[j]][3])|Km.st[[j]][3] == 1){
            ## Km[[j]] <- c(dim(Zb.st[,,ej[[j]]==1]), 1)
            dim(MU.st[[j]]) <- c(Km[[j]][1], Km.st[[j]][2])
            ZY[[j]] <- array(Z[,,ej[[j]]==1], dim = c(Km.st[[j]][1], Km.st[[j]][2], 1))
            Zm[[j]] <- array(Zb.st[,,ej[[j]]==1], dim = c(Km.st[[j]][1], Km.st[[j]][2], 1))
            ZU[[j]] <- array(Z[,,ej[[j]]==1] - MU.st[[j]], dim = c(Km.st[[j]][1], Km.st[[j]][2], sum(ej[[j]]==1)))
            ## Ym[[j]] <- array(Y[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], 1))
        } else{                ## Km.st[[j]] <- dim(Zb[,,ej[[j]]==1])
            ZY[[j]] <- Z[,,ej[[j]]==1]
            Zm[[j]] <- Zb.st[,,ej[[j]]==1]
            ZU[[j]] <- Z[,,ej[[j]]==1] - MU.st[[j]]
            ## Ym[[j]] <- Y[,,ej[[j]]==1]
        }
        ## UTA array: TRUE for upper triangle
        UTA[[j]] <- Zm[[j]]*NA
        for(k in 1:Km.st[[j]][3]) {
            UTA[[j]][,,k] <-  upper.tri(Zm[[j]][,,1])
        } 
        UTA[[j]] <- (UTA[[j]]==1)
    }
    
    
##########################################################################
    ## Likelihood computation
##########################################################################
    
    ## If sequential = FALSE compute marginal likelihood
    loglike.upper.t <- sapply(1:K[3], function(t){sum(dnorm(c(Zb.st[,,t][UTAsingle]),
                                                            mean = c(MU.state.st[[s.st[t]]][,,t][UTAsingle]),
                                                            sd=sqrt(s2.st[[s.st[t]]]), log=TRUE))})
    loglike.upper <- sum(loglike.upper.t)
   
    cat("    loglike: ", as.numeric(loglike.upper), "\n")
    
    
##########################################################################
    ## Waic
##########################################################################
    if(Waic){
        for(j in 1:ns){
            ## Waic computation
            Z.loglike.mat <- matrix(Z.loglike.array, nrow = nstore, ncol=N.upper.tri*K[3])
            ##
            ## 
            ## Z.loglike.mat[!is.finite(Z.loglike.mat)] <- 0
            Waic.out <- waic(Z.loglike.mat)$total
            rm(Z.loglike.mat)
            
            ## cat("--------------------------------------------------------------------------------------------",'\n')
            cat("    Waic: ", Waic.out[1], "\n")
            ## cat("lpd: ", Waic.out[3], "\n")
            ## cat("p_Waic: ", Waic.out[4], "\n")
            cat("----------------------------------------------",'\n')
        }
    }
    
   
    output <- lapply(MU.record, function(x){x/nss}) ## MU.record/nss
    names(output) <- "MU"
    attr(output, "title") <- "NetworkChange Posterior Sample"
    attr(output, "Z") <- Z
    attr(output, "m") <- m
    attr(output, "R") <- R
    attr(output, "U") <- U
    attr(output, "V") <- V
    attr(output, "SOS") <- SOS
    attr(output, "Umat") <- Umat
    attr(output, "Vmat") <- Vmat
    if(constant){attr(output, "bmat") <- bhat.mat}
    attr(output, "s2mat") <- s2mat
    attr(output, "Smat") <- Smat ## state output
    attr(output, "mcmc") <- nstore
    attr(output, "Waic.out") <- Waic.out
    attr(output, "prob.state") <- ps.store/nss
    ## attr(output, "loglike") <- loglike
    attr(output, "loglike") <- loglike.upper
    ## cat("elapsed time for m = ", m, " is ", proc.time() - ptm, "\n")
    return(output)
}

