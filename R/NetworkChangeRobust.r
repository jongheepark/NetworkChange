#' Changepoint analysis of a degree-corrected multilinear tensor model with t-distributed error
#'
#' NetworkChangeRobust implements Bayesian multiple changepoint models to network time series data
#' using a degree-corrected multilinear tensor decomposition method with t-distributed error

#' @param Y Reponse tensor 
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
#' @param degree.normal	A null model for degree correction. Users can choose "NULL", "eigen" or "Lsym."
#' "NULL" is no degree correction. "eigen" is a principal eigen-matrix consisting of
#' the first eigenvalue and the corresponding eigenvector. "
#' Lsym" is a modularity matrix. Default is "eigen."
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
#' 
#' @param b0 The prior mean of \eqn{\beta}. This must be a scalar. The default value is 0.
#' @param B0 The prior variance of \eqn{\beta}. This must be a scalar.  The default value is 1.
#' @param c0 = 0.1 The shape parameter of inverse gamma prior for \eqn{\sigma^2}. 
#' @param d0 = 0.1 The rate parameter of inverse gamma prior for \eqn{\sigma^2}. 
#' @param n0 = 0.1 The shape parameter of inverse gamma prior for \eqn{\gamma} of Student-t distribution. 
#' @param m0 = 0.1 The rate parameter of inverse gamma prior for \eqn{\gamma} of Student-t distribution.  
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

#' @return An mcmc object that contains the posterior sample. This object can
#' be summarized by functions provided by the coda package. The object
#' contains an attribute \code{Waic.out} that contains results of WAIC and the log-marginal
#' likelihood of the model (\code{logmarglike}). The object
#' also contains an attribute \code{prob.state} storage matrix that contains the
#' probability of \eqn{state_i} for each period
#'
#' @seealso \code{\link{NetworkStatic}}
#'
#' @references    Jong Hee Park and Yunkyun Sohn. 2020. "Detecting Structural Change
#' in Longitudinal Network Data." \emph{Bayesian Analysis}. Vol.15, No.1, pp.133-157.
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
#'    G <- 100 ## only 100 mcmc scans to save time
#'    ## Fit models
#'    out1 <- NetworkChangeRobust(Y, R=2, m=1, mcmc=G, burnin=G, verbose=G)
#'    ## plot latent node positions
#'    plotU(out1)
#'    ## plot layer-specific network generation rules
#'    plotV(out1)
#'    }

##################################################################################
NetworkChangeRobust <- function(Y, R=2, m=1, initial.s = NULL,  
                                mcmc=100, burnin=100, verbose=0, thin  = 1,                         
                                degree.normal="eigen", 
                                UL.Normal = "Orthonormal",
                                plotUU = FALSE, plotZ = FALSE,
                                b0 = 0, B0 = 1, c0 = NULL, d0 = NULL, n0 = 2, m0 = 2,
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
    Z <- Y
    nss <- 0
    K <- dim(Y)  
    Time <- K[3]
    P  <-  trans.mat.prior(m=m, n=Time, a = 0.9, b= 0.1)
    A0  <-  trans.mat.prior(m=m, n=Time, a = a, b = b)
    
    ## if (is.null(initial.V)){
    out <- startUV(Z, R, K)
    initial.U <- out[[1]]
    V <- out[[2]]
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
     
    nodenames <- dimnames(Y)[[1]]
    
    ## unique values of Y
    uy <- sort(unique(c(Y)))
    
    ## degree normalization
    if(degree.normal == "eigen"){
        ## all the time
        for(k in 1:K[3]){
            ee <- eigen(Y[,,k])
            Z[,,k] <- Y[,,k] - ee$values[1] * outer(ee$vectors[,1], ee$vectors[,1])
            diag(Z[,,k]) <-  0
        }
    }
    
    ## if Modularity
    if(degree.normal == "Modul"){
        gamma.par = 1
        for(k in 1:K[3]){
            Yk <- as.matrix(Y[,,k])
            yk <- as.vector(apply(Yk, 2, sum))
            ym <- sum(yk)
            Z[,,k] <- Yk - gamma.par*(yk%o%yk)/ym
            diag(Z[,,k]) <-  0
        }
    }
    
    
    ## initialize beta
    bhat <- mean(c(Z))
    X <- array(1, dim=c(K, 1))
    p <- dim(X)[4]
    XtX <- prod(K) ## matrix(sum(X^2), p, p)
    rm(X)
    Zb <- Z - bhat
    
    
    ## eigen decomposition
    ## VM is time specific eigenvalues
    
    if (is.null(initial.s)){
        s <- startS(Z, Time, m, initial.U, V, s2=1, R)
    } else{
        s <- initial.s
    }
    
    ## holder 
    ## Zm is a state-specific holder of ZE = Z - bhat
    ## Zm[[1]] is a subset of Z pertaining to state 1
    ## ZU = Z - ULU
    ## ZY is original Z separated by state
    UTA <- Km <- Zm <- ZY <- ZU <- ej <- U <- MU <- MU.state <- Xm <- Vm <- as.list(rep(NA, ns))
    EEt <- tEE <- lambda <- n2 <- ZEE <- as.list(rep(NA, ns))
    ps.store <- matrix(0, Time, ns)
    
    ## given the state vector, initialize regime specific U and Vm
    for (j in 1:ns){
        ej[[j]] <- as.numeric(s==j)
        Zm[[j]] <- Zb[,,ej[[j]]==1] 
        tmp <- eigen(apply(Zm[[j]], c(1,2), mean))
        d2m <- abs(tmp$val)
        U0 <- tmp$vec[, order(d2m, decreasing=TRUE) ]
        U[[j]] <- matrix(U0[, 1:R], nrow=nrow(U0), ncol=R)
        Vm[[j]] <- matrix(V[ej[[j]] == 1, ], sum(s==j), R)
        ej[[j]] <- as.numeric(s==j)
        Km[[j]] <- dim(Z[,,ej[[j]]==1])
        if(is.na(Km[[j]][3])) Km[[j]][3] <- 1
     }
    ## V <- Reduce(rbind, Vm)
    
    ## initialize MU and MU.state
    ## MU is regime-specific mean matrix, the length of which depends on regime length
    ## MU.state is a full-length mean matrix for state sampling
    for (j in 1:ns){
        MU[[j]] <-  M.U(list(U[[j]],U[[j]], Vm[[j]]))
        MU.state[[j]] <-  M.U(list(U[[j]],U[[j]],V))
    }
    
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
        lambda[[j]] <- rep(1, sum(s==j))
    }
    bhat.mat <- rep(NA, nstore)
    Vmat <- matrix(NA, nstore, R*K[3])
    Smat <- matrix(NA, nstore, K[3])
    
    ## loglike holder
    N.upper.tri <- K[1]*(K[1]-1)/2
    ## Z.loglike <- matrix(NA, mcmc, K[3])
    ## Z.loglike <- as(matrix(NA, mcmc, K[3]), "mpfr")

    Zt <- matrix(NA,  Time,  N.upper.tri)
    UTAsingle <-  upper.tri(Z[,,1])          
    Waic.out <- NA
    SOS <- 0
    
    ## lambda initialization
    n1 <- n0 + 1
    for(j in 1:ns){
        n2[[j]] <- m0
        lambda[[j]] <- rgamma(Km[[j]][3], n1/2, n2[[j]]/2)
    }

    
    if(verbose !=0){
        cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n")
        cat("\t NetworkChangeRobust Sampler Starts! \n")
        ## cat("\t function called: ")
        ## print(call)
        cat("\t degree normalization: ", degree.normal, "\n")
        cat("\t initial states: ", table(s), "\n")
        cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n")
    }
    SOS.random = TRUE
#############################################################
    ## MCMC loop starts!
#############################################################
    for(iter in 1:totiter) {
        if(iter > burnin){
            random.perturb = FALSE
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
                ## Xm[[j]] <- array(X[,,ej[[j]]==1, ], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1), p))
                ## ZU[[j]] <- array(Z[,,ej[[j]]==1] - MU[[j]], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1)))
                 ## Ym[[j]] <- array(Y[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], 1))
            } else{
                ZY[[j]] <- Z[,,ej[[j]]==1]
                Zm[[j]] <- Zb[,,ej[[j]]==1]
                ## Xm[[j]] <- array(X[,,ej[[j]]==1, ], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1), p))
                ## ZU[[j]] <- Z[,,ej[[j]]==1] - MU[[j]]
                ## Ym[[j]] <- Y[,,ej[[j]]==1]
            }
           
            ## return the right dimension info
            Km[[j]] <- dim(Zm[[j]])
            if(is.na(Km[[j]][3])) {
                Km[[j]][3] <- 1
                Zm[[j]] <- array(Zm[[j]], dim = c(dim(Zm[[j]])[1], dim(Zm[[j]])[2], 1))
            }
            ## UTA array: TRUE for upper triangle
            UTA[[j]] <- Zm[[j]]*NA
            for(k in 1:Km[[j]][3]) {
                UTA[[j]][,,k] <-  upper.tri(Zm[[j]][,,1])
            } 
            UTA[[j]] <- (UTA[[j]]==1)

            ## lambda reset
            lambda[[j]] <- united.lambda[ej[[j]]==1]
        }

       
        ## Step 2. update U
        ## cat("\n---------------------------------------------- \n ")
        ## cat("Step 2. update U \n")
        ## cat("\n---------------------------------------------- \n ")
        ## updateUm <- function(ns, U, V, R, Zm, Km, ej, s2, eU, iVU, UL.Normal){
        for(j in 1:ns){
            Vj <-  matrix(V[ej[[j]] == 1, ], nrow=sum(ej[[j]]), ncol=R)
            ## cat("s = ", table(s), "\n")
            ## cat("Vj = ", dim(Vj), "\n")
            ## cat("lambda[[j]] = ", length(lambda[[j]]), "\n")
            for(i in sample(Km[[j]][1])){
                Ui <- U[[j]]
                Ui[i,] <- 0
                VU <-  aperm(array(apply(Ui,1,"*",t(Vj*lambda[[j]])), dim=c(R, Km[[j]][3], Km[[j]][1])), c(3,2,1))
                zi <- Zm[[j]][i,,]
                L <-  apply(VU*array(rep(zi,R), dim=c(Km[[j]][1], Km[[j]][3], R)), 3, sum) 
                Q <-  (t(Ui)%*%Ui) * (t(Vj*lambda[[j]])%*%Vj)
                cV <- solve(Q/s2[j] + iVU[[j]] )
                cE <- cV%*%( L/s2[j] + iVU[[j]]%*%eU[[j]])
                U[[j]][i,] <- rMVNorm(1, cE, cV ) 
            }
        }
        ## UL normalization
        if (UL.Normal == "Normal"){
            for(j in 1:ns){
                U[[j]] <- Unormal(U[[j]])
            }
        }else if(UL.Normal == "Orthonormal"){
            for(j in 1:ns){
                U[[j]] <- GramSchmidt(U[[j]])
            }
        }
        ## return(U)
        ## }
        
        ## U <- updateUm(ns, U, V, R, Zm, Km, ej, s2, eU, iVU, UL.Normal)
        
        ## Step 3. update V
        ## cat("\n---------------------------------------------- \n ")
        ## cat("Step 3. update V \n")
        ## cat("\n---------------------------------------------- \n ")
        ## Vm <- updateVm(ns, U, V, Zm, Km, R, s2, eV, iVV, UTA)
        Vm <- as.list(rep(NA, ns))
        for(j in 1:ns){
            Uj <- U[[j]]
            Zj <- Zm[[j]]
            Zj[!UTA[[j]]] <- 0            
            Q <- UU <- ZEP <- L <- cV <- cE <- NA
            if(R == 1){
                Q <- ((t(Uj)%*%Uj)^2 - matrix(sum(Uj^4), R, R))/2
            } else{
                Q <- ((t(Uj)%*%Uj)^2 -
                                   matrix(apply(apply(Uj^2,1,function(x){x%*%t(x)}), 1, sum), R, R))/2
            }
            UU <- aperm(array(apply(Uj,1,"*",t(Uj)),
                              dim=c(R, Km[[j]][1], Km[[j]][1]) ),c(2,3,1))
            ZEP <- aperm(Zj, c(3,1,2))
            ZUU <- array(apply(UU,3,function(x){apply(ZEP,1,"*",x)}),
                         dim=c(Km[[j]][1], Km[[j]][1], Km[[j]][3], R))
            L <- apply(ZUU, c(3,4),sum)
            ## cV <- solve(Q/s2[j] + iVV[[j]])
            ## weight cV by the sum of lambda at state j
            cV <- solve(Q/s2[j] + iVV[[j]])
            cE <- (L/s2[j] + rep(1, Km[[j]][3])%*%t(eV[[j]])%*%iVV[[j]])%*%cV    
            Vm[[j]] <-  rmn(cE, diag(Km[[j]][3]), cV)
        }
        V <- Reduce(rbind, Vm)

        ## update MU
        for(j in 1:ns){
            ## MU is shorter than MU.state. MU.state is a full length.
            MU[[j]] <- M.U(list(U[[j]],U[[j]],Vm[[j]]))
            MU.state[[j]] <- M.U(list(U[[j]],U[[j]],V))
            ZU[[j]] <- ZY[[j]] - MU[[j]]
        }
        
        ## Step 4. update lambda
        n1 <- n0 + 1        
        for(j in 1:ns){
            ## cat("Zm[[j]] = ", dim(Zm[[j]]), "\n")
            ## cat("MU[[j]] = ", dim(MU[[j]]), "\n")
            ZEE[[j]] <- Zm[[j]] - MU[[j]]
            tEE[[j]] <- apply(ZEE[[j]], 3, c) ## NN by T1
            EEt[[j]] <- sapply(1:Km[[j]][3], function(gg){sum(tEE[[j]][,gg]^2)/s2[[j]]}) 
            n2[[j]] <- m0 + EEt[[j]]
            lambda[[j]] <- rgamma(Km[[j]][3], n1/2, n2[[j]]/2)
        }

        
        ## Step 5. update s2
        ## cat("\n---------------------------------------------- \n ")
        ## cat("Step 4. update s2 \n")
        ## cat("\n---------------------------------------------- \n ")
        ## s2 <- updates2m(ns, Zm, MU, c0, d0, Km)
        for(j in 1:ns){
            ## ZEE[[j]] <- Zm[[j]] - MU[[j]]
            ## EE <- c(ZEE[[j]])                    
            ## tEE <- apply(ZEE[[j]], 3, c) ## NN by T1
            s2[j] <- 1/rgamma(1, (c0+prod(Km[[j]]))/2, (d0+ sum(EEt[[j]]*lambda[[j]]))/2)
        }  
        
        
        ## update constant bhat
        bhat <- updatebm(ns, K, s, s2, B0, p, ZU)
        Zb <- Z - bhat
      
        ## update hierarchical parameters
        ## hierarchical parameters for U
        for(j in 1:ns){
            SS <-  t(U[[j]]) %*% U[[j]]## (Km[[j]][1]-1)*cov(U[[j]]) + Km[[j]][1]*msi/(Km[[j]][1]+1)
            for(r in 1:R){
                iVU[[j]][r,r] <- 1/rgamma(1, (u0 + K[1])/2, (u1+ SS[r,r])/2)
            }
            eU[[j]] <- c(rMVNorm(1,apply(U[[j]],2,sum)/(Km[[j]][1]+1), solve(iVU[[j]])/(Km[[j]][1]+1)))
        }
        
        ## hierarchical parameters for V
        ## V for state j only
        for(j in 1:ns){
            Vs <- matrix(Vm[[j]], nrow=sum(ej[[j]]), ncol=R)
            SS <-  t(Vs)%*%Vs
            for(r in 1:R){
                iVV[[j]][r,r] <- 1/rgamma(1, (v0 + Km[[j]][3])/2, (v1 + SS[r,r])/2)
            }
            eV[[j]] <- c(rMVNorm(1,apply(Vs, 2, sum)/(Km[[j]][3]+1),
                                                 solve(iVV[[j]])/(Km[[j]][3]+1)))      
        }
        
        

       ## Zb = Z - bhat
        ## Zm is regime specific Zb
        ## ZY = Z
        ## ZU = ZY - MU
        
        ## Step 1. update s
        ## cat("\n---------------------------------------------- \n ")
        ## cat("Step 5. update s \n")
        ## cat("\n---------------------------------------------- \n ")
        ## state.out <- updateS(iter, s, V, m, Zb, Zt, Time, MU.state, P, s2,
        ##                     N.upper.tri, random.perturb)
        ## state.out <- updateS(iter, s, V, m, Zb, Zt, Time, fast,
        ##                       MU.state, P, s2, local.type, logistic.tune, N.upper.tri, sticky)
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
        ## state.out <- ULUstateSample(m=m, s=s, ZMUt=ZMUt, s2=s2, P=P, random.perturb)
        united.lambda <- unlist(lambda)
        density.log <- as(matrix(unlist(sapply(1:T, function(t){
            lapply(ZMUt, function(x){sum(dnorm(x[t,], 0,
                                               sd = sqrt(s2[s[t]]/united.lambda[t]), log=TRUE))})})), ns, Time), "mpfr")
        Fmat   <-  as(matrix(NA, Time, m+1), "mpfr")     # storage for the Filtered probabilities
        pr1 <-  as(c(1,rep(0, m)), "mpfr")         # initial probability Pr(s=k|Y0, lambda)
        unnorm.pstyt <- py <- as(rep(NA, m+1), "mpfr")
        for (t in 1:Time){
            if(t==1) {
                pstyt1 = pr1
            }else {
                pstyt1 <- Fmat[t-1,]%*%P
            }                
            ## par(mfrow=c(1,2));
            unnorm.pstyt    <- pstyt1*exp(density.log[,t])       
            Fmat[t,]  <-  unnorm.pstyt/sum(unnorm.pstyt) # Pr(st|Yt)
            ## cat("\nunnorm.pstyn at t = ", t, " is ", unnorm.psty, "\n")
        }
        Fmat <- matrix(as.numeric(Fmat), nrow=Time, m+1)
        s      <-  matrix(1, Time, 1)   ## holder for state variables
        ps     <-  matrix(NA, Time, m+1) ## holder for state probabilities
        ps[Time,] <-  Fmat[Time,]              ## we know last elements of ps and s
        s[Time,1] <-  m+1
        ## t      <-  T-1
        for(t in (Time-1):2){
            ## while (t>=1){
            st     <-  s[t+1]
            unnorm.pstyn   <-  Fmat[t,]*P[,st]
            ## cat("\nunnorm.pstyn at t = ", t, " is ", unnorm.pstyn, "\n")
            if(sum(unnorm.pstyn) == 0){
                ## if unnorm.pstyn is all zero, what to do?
                cat("F", Fmat[t,]," and P", P[,st]," do not match at t = ", t, "\n")
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
        new.SOS <- FALSE
        if(SOS.random & sum(table(s) == 1)){
            s <- sort(sample(1:ns, T, replace=TRUE, prob=rep(1/ns, ns))) ## 
            ## cat("\n A single observation state is sampled and the latent state sampled randomly.\n")
            if(length(unique(s)) != ns){
                s <- sort(rep(1:ns, length=T))
            }
            new.SOS <- TRUE
        }
        ## s <- state.out$s
        ## ps <- state.out$ps
        SOS <- SOS + new.SOS
        
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
        P <- updateP(s, ns, P, A0)

        
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
            ## s2mat[iter-burnin] <- s2
            for(j in 1:ns){
                MU.record[[j]] <- MU.record[[j]] + MU.state[[j]]
                s2mat[[j]][(iter-burnin)/thin] <- s2[j]
                ## bhat.mat[[j]][iter-burnin] <- bhat[[j]]
                ## likelihood for upper triangle array
                Umat[[j]][(iter-burnin)/thin, ] <- as.vector(U[[j]])
                eUmat[[j]][(iter-burnin)/thin, ] <- as.vector(eU[[j]])
                iVUmat[[j]][(iter-burnin)/thin, ] <- as.vector(iVU[[j]])
                eVmat[[j]][(iter-burnin)/thin, ] <- as.vector(eV[[j]])
                iVVmat[[j]][(iter-burnin)/thin, ] <- as.vector(iVV[[j]])
            }
            Vmat[(iter-burnin)/thin, ] <- as.vector(V)
            bhat.mat[(iter-burnin)/thin] <- bhat
            Smat[(iter-burnin)/thin, ] <- s
            Pmat[(iter-burnin)/thin, ] <- diag(P)
            ps.store <- ps.store + ps
            ## if(DIC|Waic){
            ##     d <- sapply(1:K[3], function(t){dnorm(c(Zb[,,t][UTAsingle]),
            ##                                           mean = c(MU.state[[s[t]]][,,t][UTAsingle]),
            ##                                           sd=sqrt(s2[[s[t]]]), log=TRUE)})
            ##     Z.loglike.array[(iter-burnin)/thin, ,] <- d
            ## }
        }
        
    }## end of MCMC loop
    ## sort( sapply(ls(),function(x){object.size(get(x))}))
 
    #########################################################################
    ## flip V and U based on the size of mse(V)
    ## so that the first principal axis comes at the first column of V and U
    ## In that way, the interpretation of cp results makes sense.
    #########################################################################
    
    output <- lapply(MU.record, function(x){x/nss}) ## MU.record/nss
    names(output) <- "MU"
    attr(output, "title") <- "NetworkChangeRobust Posterior Sample"
    attr(output, "Z") <- Z
    attr(output, "m") <- m
    attr(output, "R") <- R
    attr(output, "U") <- U
    attr(output, "V") <- V
    attr(output, "SOS") <- SOS
    attr(output, "Umat") <- Umat
    attr(output, "Vmat") <- Vmat
    attr(output, "bmat") <- bhat.mat
    attr(output, "s2mat") <- s2mat
    attr(output, "Smat") <- Smat ## state output
    attr(output, "mcmc") <- nstore
    ## attr(output, "DIC") <- sum(as.numeric(Z.DIC))
    ## attr(output, "Waic.out") <- Waic.out
    attr(output, "prob.state") <- ps.store/nss
    ## attr(output, "loglike") <- loglike
    ## attr(output, "loglike") <- loglike.upper
    ## attr(output, "logmarglike") <- logmarglike
    ## attr(output, "logmarglike") <- logmarglike.upper 
    ## cat("elapsed time for m = ", m, " is ", proc.time() - ptm, "\n")
    return(output)
}

