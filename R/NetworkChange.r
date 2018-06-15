#' Changepoint analysis of a degree-corrected multilinear tensor model
#'
#' NetworkChange implements Bayesian multiple changepoint models to network time series data
#' using a degree-corrected multilinear tensor decomposition method

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
#' @param constant If \code{constant = TRUE}, constant parameter is sampled
#' and saved in the output as attribute \code{bmat}. Default is \code{constant = FALSE}.
#' 
#' @param b0 The prior mean of \eqn{\beta}. This must be a scalar. The default value is 0.
#' @param B0 The prior variance of \eqn{\beta}. This must be a scalar.  The default value is 1.
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

#' @param marginal If \code{marginal = TRUE}, the log marignal likelihood is computed using the method of Chib (1995).
    
#' @param DIC If \code{DIC = TRUE}, the deviation information criterion is computed.
    
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
#'    out0 <- NetworkStatic(Y, R=2, mcmc=G, burnin=G, verbose=G, Waic=TRUE)
#'    out1 <- NetworkChange(Y, R=2, m=1, mcmc=G, burnin=G, verbose=G, Waic=TRUE)
#'    out2 <- NetworkChange(Y, R=2, m=2, mcmc=G, burnin=G, verbose=G, Waic=TRUE)
#'    out3 <- NetworkChange(Y, R=2, m=3, mcmc=G, burnin=G, verbose=G, Waic=TRUE)
#'    outlist <- list(out0, out1, out2, out3)
#'
#'    ## The most probable model given break number 0 to 3 and data is out1 according to WAIC 
#'    WaicCompare(outlist)
#'
#'    ## plot latent node positions
#'    plotU(out1)
#'  
#'    ## plot layer-specific network generation rules
#'    plotV(out1)
#'    }

##################################################################################
NetworkChange <- function(Y, R=2, m=1, initial.s = NULL,  
                          mcmc=100, burnin=100, verbose=0, thin  = 1,                         
                          degree.normal="eigen", 
                          UL.Normal = "Orthonormal",
                          DIC = FALSE, Waic=FALSE, marginal = FALSE, 
                          plotUU = FALSE, plotZ = FALSE, constant = FALSE,
                          b0 = 0, B0 = 1, c0 = NULL, d0 = NULL,
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

    if(R==1 & UL.Normal == "Orthonormal"|| R==1 & UL.Normal == "Normal"){
        stop("If R=1, please set UL.Normal=FALSE.")
    }

    
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
    if(constant){
        bhat <- mean(c(Z))
        Zb <- Z - bhat
    }else{
        bhat = 0
        Zb <- Z
    }
    X <- array(1, dim=c(K, 1))
    p <- dim(X)[4]
    XtX <- prod(K) ## matrix(sum(X^2), p, p)
    rm(X)
    
    
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
    }
    if(constant){
        bhat.mat <- rep(NA, nstore)
    }
    Vmat <- matrix(NA, nstore, R*K[3])
    Smat <- matrix(NA, nstore, K[3])
    
    ## loglike holder
    N.upper.tri <- K[1]*(K[1]-1)/2
    ## Z.loglike <- matrix(NA, mcmc, K[3])
    ## Z.loglike <- as(matrix(NA, mcmc, K[3]), "mpfr")
    if(DIC|Waic){
        Z.loglike.array <- array(NA, dim=c(nstore, N.upper.tri, K[3]))
    }
    logmarglike <- loglike <- logmarglike.upper <- loglike.upper <- NA

    Zt <- matrix(NA,  Time,  N.upper.tri)
    UTAsingle <-  upper.tri(Z[,,1])          
    Waic.out <- NA
    SOS <- 0
    
    if(verbose !=0){
        cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n")
        cat("\t NetworkChange Sampler Starts! \n")
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
        if(iter > burnin){
            random.perturb = FALSE
        }
        ## Zb = Z - bhat
        ## Zm is regime specific Zb
        ## ZY = Z
        ## ZU = ZY - MU
        
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
        U <- updateUm(ns, U, V, R, Zm, Km, ej, s2, eU, iVU, UL.Normal)
        
        ## Step 3. update V
        ## cat("\n---------------------------------------------- \n ")
        ## cat("Step 3. update V \n")
        ## cat("\n---------------------------------------------- \n ")
        Vm <- updateVm(ns, U, V, Zm, Km, R, s2, eV, iVV, UTA)
        V <- Reduce(rbind, Vm)

        ## update MU
        for(j in 1:ns){
            ## MU is shorter than MU.state. MU.state is a full length.
            MU[[j]] <- M.U(list(U[[j]],U[[j]],Vm[[j]]))
            MU.state[[j]] <- M.U(list(U[[j]],U[[j]],V))
            ZU[[j]] <- ZY[[j]] - MU[[j]]
        }
        
        ## Step 4. update s2
        ## cat("\n---------------------------------------------- \n ")
        ## cat("Step 4. update s2 \n")
        ## cat("\n---------------------------------------------- \n ")
        s2 <- updates2m(ns, Zm, MU, c0, d0, Km)
       
        ## update constant bhat
        if(constant){
            bhat <- updatebm(ns, K, s, s2, B0, p, ZU)
            Zb <- Z - bhat
        }
        
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
        
        
        ## Step 5. update s
        ## cat("\n---------------------------------------------- \n ")
        ## cat("Step 5. update s \n")
        ## cat("\n---------------------------------------------- \n ")
        state.out <- updateS(iter, s, V, m, Zb, Zt, Time, MU.state, P, s2,
                             N.upper.tri, random.perturb)
        ## state.out <- updateS(iter, s, V, m, Zb, Zt, Time, fast,
        ##                       MU.state, P, s2, local.type, logistic.tune, N.upper.tri, sticky)
        s <- state.out$s
        ps <- state.out$ps
        SOS <- SOS + state.out$SOS
 
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
            if(constant){ bhat.mat[(iter-burnin)/thin] <- bhat}
            Smat[(iter-burnin)/thin, ] <- s
            Pmat[(iter-burnin)/thin, ] <- diag(P)
            ps.store <- ps.store + ps
            if(DIC|Waic){
                d <- sapply(1:K[3], function(t){dnorm(c(Zb[,,t][UTAsingle]),
                                                      mean = c(MU.state[[s[t]]][,,t][UTAsingle]),
                                                      sd=sqrt(s2[[s[t]]]), log=TRUE)})
                Z.loglike.array[(iter-burnin)/thin, ,] <- d
            }
        }
        
    }## end of MCMC loop
    ## sort( sapply(ls(),function(x){object.size(get(x))}))
    ## rm(list = c('ZEP','Zj', 'UTA', "MUt", "ZMUt", "Y", "Zm", "ZUU", "ZU", "ZY", "MU.state"))
    ##               ZEP                 Zj                UTA                MUt 
    ##          ZMUt                  Y                  Z                 Zb 
    ##           ZUU                 MU                ZEE                 Zm 
    ##            ZU                 ZY               MU.record           MU.state 



    
##########################################################################
    ## Model Diagnostics
##########################################################################
    ## DIC computation = 2*D_average - D_theta^hat
    ## D_theta^hat: posterior estimates
    s.st <- ceiling(apply(Smat, 2, median))
    ## MU.record.st <-  lapply(MU.record, function(x){x/nss})## MU.record[[j]]/nss      
    M.st <-  lapply(MU.record, function(x){x/nss})## MU.record[[j]]/nss      
    if(constant){
        bhat.st <- mean(bhat.mat) ## mean((trueY - M.st)^2 )
    }else{
        bhat.st <- 0
    }
        
    P.st <- apply(Pmat, 2, mean)
    Km.st <- Km
    for(j in 1:ns){
        Km.st[[j]][3] <- sum(s.st == j)
    }
    if(constant){
        Zb.st <- Z - bhat.st
    }else{
        Zb.st <- Z
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
    if(sequential){
        loglike.T <- loglike.upper.T <- matrix(NA, K[3], ns)
        for(j in 1:ns){
            loglike.upper.T[,j] <- sapply(1:K[3], function(t){sum(dnorm(c(Zb.st[,,t][UTAsingle]),
                                                                        mean = c(MU.state.st[[j]][,,t][UTAsingle]),
                                                                        sd=sqrt(s2.st[[j]]), log=TRUE))})
            ## loglike.T[,j] <- sapply(1:K[3], function(t){sum(dnorm(c(Zb.st[,,t]),
            ##                                                      mean = c(MU.state.st[[j]][,,t]),
            ##                                                      sd=sqrt(s2.st[[j]]), log=TRUE))})
        }
        ## sequential estiamte of loglike does not make a difference
        MUt <- list()
        for (t in 1:Time){
            Zt[t, ] <- c(Zb.st[, , t][upper.tri(Zb.st[, , t])])
        }
        for(j in 1:ns){
            MUt[[j]] <- matrix(NA,  Time, N.upper.tri)
            for (t in 1:Time){
                MUt[[j]][t, ] <- c((MU.state.st[[j]][, , t])[upper.tri(MU.state.st[[j]][, , t])])
            }
        }
        ZMUt <- as.list(rep(NA, ns))
        for(j in 1:(m+1)){
            ZMUt[[j]] <- Zt - MUt[[j]]
        }
        
        density.log <- as(matrix(unlist(sapply(1:Time, function(t){
            lapply(ZMUt, function(x){sum(dnorm(x[t,], 0, sd = sqrt(s2.st), log=TRUE))})})), ns, Time), "mpfr")
        F   <-  as(matrix(NA, Time, m+1), "mpfr")     # storage for the Filtered probabilities
        pr1 <-  as(c(1,rep(0, m)), "mpfr")         # initial probability Pr(s=k|Y0, lambda)
            unnorm.pstyt <- py <- as(rep(NA, m+1), "mpfr")
        pstyt1 <-  as(matrix(NA, Time, m+1), "mpfr")
        for (t in 1:Time){
            if(t==1) {
                pstyt1[t, ] = pr1
            }else {
                pstyt1[t, ] <- F[t-1,]%*%P
            }                
            ## par(mfrow=c(1,2));
            unnorm.pstyt    <- pstyt1[t, ]*exp(density.log[,t])       
            F[t,]  <-  unnorm.pstyt/sum(unnorm.pstyt) # Pr(st|Yt)
        }
        pst <- matrix(as.numeric(pstyt1), nrow=Time, m+1)
        ## loglike.t <-  apply(loglike.T*(pst), 1, sum)
        loglike.upper.t <- apply(loglike.upper.T*(pst), 1, sum)
        ## loglike <- sum(loglike.t)
        loglike.upper <- sum(loglike.upper.t)
    }else{
        ## If sequential = FALSE compute marginal likelihood
        loglike.upper.t <- sapply(1:K[3], function(t){sum(dnorm(c(Zb.st[,,t][UTAsingle]),
                                                                mean = c(MU.state.st[[s.st[t]]][,,t][UTAsingle]),
                                                                sd=sqrt(s2.st[[s.st[t]]]), log=TRUE))})
        ## loglike.t <- sapply(1:K[3], function(t){sum(dnorm(c(Zb.st[,,t]),
        ##                                                  mean = c(MU.state.st[[s.st[t]]][,,t]),
        ##                                                  sd=sqrt(s2.st[[s.st[t]]]), log=TRUE))})
        ## loglike <- sum(loglike.t)
        loglike.upper <- sum(loglike.upper.t)
    }
    
    ## cat("    loglike: ", as.numeric(loglike), "\n")
    if(verbose !=0){
        cat("    loglike : ", as.numeric(loglike.upper), "\n")
    }
 
    
##########################################################################
    ## DIC and Waic
##########################################################################
    logpy <- pDIC <- Z.DIC <- ## pDIC.gelman.alt <- Z.DIC.gelman.alt <-
        sumColMean <- sumColVar <- rep(NA, ns)
    
    
    if(DIC|Waic){
        for(j in 1:ns){
            ## DIC3 by Celeux, Forbes, Robert, and Titterington p.655
            ## is equivalent to Gelman, Hwang, Vehtari's pDIC1
            ## logpy <- sum(log(dnorm(c(ZU), mean = c(MUU), sqrt(s2.st))), na.rm=TRUE)
            ## pDIC <- 2*(logpy - sum(colMeans(Z.loglike.mat)))
            ## Z.DIC <- -2*logpy + 2*pDIC
            if(length(which(s.st == j)) == 1){
                t <- which(s.st == j)
                logpy[j] <- sum(dnorm(c(Zb.st[,,t][UTAsingle]),
                                      mean = c(MU.state.st[[s.st[t]]][,,t][UTAsingle]),
                                      sd=sqrt(s2.st[[s.st[t]]]), log=TRUE))
            }else{
                loglike.upper.t <- sum(sapply(which(s.st == j), function(t){sum(dnorm(c(Zb.st[,,t][UTAsingle]),
                                                        mean = c(MU.state.st[[s.st[t]]][,,t][UTAsingle]),
                                                        sd=sqrt(s2.st[[s.st[t]]]), log=TRUE))}))
                logpy[j] <- sum(loglike.upper.t)
            }
            Z.loglike.mat <- Z.loglike.array[,,which(s.st == j)]
            Z.loglike.mat[!is.finite(Z.loglike.mat)] <- 0
            dim(Z.loglike.mat) <- c(nstore, N.upper.tri*sum(ej[[j]]))
            
            
            sumColMean[j] <- sum(colMeans(Z.loglike.mat))
            pDIC[j] <-  2*(logpy[j] - sumColMean[j])
            Z.DIC[j] <- -2*logpy[j] + 2*pDIC[j]
            
            ## sumColVar[j] <- sum(colVars(Z.loglike.mat))
            ## pDIC.gelman.alt[j] <-  sumColVar[j]
            ## Z.DIC.gelman.alt[j] <- -2*logpy[j] + 2*pDIC.gelman.alt[j]           
        }
    
        if(DIC == TRUE){    
            cat("----------------------------------------------",'\n')
            ## cat("    log p(y|theta.Bayes): ", sum(as.numeric(logpy)), "\n")
            cat("    DIC: ", sum(as.numeric(Z.DIC)), "\n")
            ## cat("    pDIC: ", sum(as.numeric(pDIC)), "\n")
            ## cat("    DIC.gelman.alt: ", sum(as.numeric(Z.DIC.gelman.alt)), "\n")
            ## cat("    pDIC.gelman.alt: ", sum(as.numeric(pDIC.gelman.alt)), "\n")
        }
        
        if(Waic == TRUE){
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
    

##########################################################################
    ## Marginal Likelihood computation
##########################################################################
    ## if (marginal == TRUE& min(table(s.st))<=1){
        ##     cat("\n\n The model is singular. Fit a model with smaller number of break. \n\n")        
        ## } else{
        ## if (marginal == TRUE & min(table(s.st))>1 ){
    if (marginal == TRUE){
        ## Hierarchical parameters must be properly considered.
        ## 1. p(eu.st|Z)
        ## 2. p(iVU.st|Z, eu.st) = p()
        ## 3. p(ed.st|Z, eu.st, iVU.st)
        ## 4. p(iVD.st|Z, eu.st, iVU.st, ed.st)
        ## 5. same
        
        ##
        ## Marginal Step 1: p(eU.st|Z)
        ##
        Zb <- Z
        density.eU.holder <- matrix(NA, reduce.mcmc)
        for (g in 1:reduce.mcmc){
            density.eU.temp.holder <- rep(NA, ns)
            for(j in 1:ns){
                ej[[j]]  <-  as.numeric(Smat[g,]==j)
                if(constant){ Zb <- Z - bhat.mat[g]}               
                if(is.na(Km[[j]][3])|Km[[j]][3] == 1){
                    Km[[j]] <- c(dim(Zb[,,ej[[j]]==1]), 1)
                    ## ZY[[j]] <- array(Z[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], 1))
                    Zm[[j]] <- array(Zb[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], 1))
                    ## Ym[[j]] <- array(Y[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], 1))
                } else{
                    Km[[j]] <- dim(Zb[,,ej[[j]]==1])
                    ## ZY[[j]] <- Z[,,ej[[j]]==1]
                    Zm[[j]] <- Zb[,,ej[[j]]==1]
                    ## Ym[[j]] <- Y[,,ej[[j]]==1]
                }
    
                ## hierarchical parameters for U
                Vj <- matrix(matrix(Vmat[g,], K[3], R)[Smat[g,] == j, ], nrow=sum(Smat[g,] == j),
                             ncol=R)
                Uj <- matrix(Umat[[j]][g,], K[1], R)
                SS <- t(Uj)%*%Uj ## (Km[[j]][1]-1)*cov(Uj) + Km[[j]][1]*msi/(Km[[j]][1]+1)
                ## invSS <- chol2inv(chol(SS))
                ## iVU[[j]] <- rwish(invSS, psi.U0 + Km[[j]][1])
                for(r in 1:R){
                    iVU[[j]][r,r] <- 1/rgamma(1, (u0 + K[1])/2, (u1 + SS[r,r])/2)
                }
                mu <- apply(Uj,2,sum)/(Km[[j]][1]+1)
                sigma <- solve(iVU[[j]])/(Km[[j]][1]+1) 
                density.eU.temp.holder[j] <- dmvnorm(eU.st[[j]], matrix(mu, 1, R), sigma, log=TRUE)
            }
            density.eU.holder[g] <- sum(density.eU.temp.holder)
        }
        density.eU <- log(mean(exp(density.eU.holder)))
        if(abs(density.eU) == Inf){
            ## cat("    Precision reinforced! \n")
            ## print(density.eU.holder)
            density.eU <- as.numeric(log(mean(exp(mpfr(density.eU.holder, precBits=53)))))
        }
        if(abs(density.eU) == Inf){
            stop("The density of a posterior ordinate reaches a computational limit and marginal likelihood computation is halted.\n")
        }
        
        cat("\n---------------------------------------------- \n ")
        cat("Marignal Likelihood Computation Step 1 \n")
        cat("    density.eU: ", as.numeric(density.eU), "\n")
        cat("---------------------------------------------- \n ")

        ##
        ## Marginal Step 2: p(iVU.st|Z, eu.st)
        ##
        density.iVU.holder <- matrix(NA, reduce.mcmc)
        for (g in 1:reduce.mcmc){
            density.iVU.temp.holder <- rep(NA, ns)
            ## Step 1. update ej, Km, Zm
            for (j in 1:ns){
                ej[[j]] <- as.numeric(s==j)
                Km[[j]] <- dim(Zb[,,ej[[j]]==1])
                if(is.na(Km[[j]][3])){
                    ZY[[j]] <- array(Z[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1)))
                    Zm[[j]] <- array(Zb[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1)))
                } else{
                    ZY[[j]] <- Z[,,ej[[j]]==1]
                    Zm[[j]] <- Zb[,,ej[[j]]==1]
                }
                Km[[j]] <- dim(Zm[[j]])
                
                ## UTA array: TRUE for upper triangle
                UTA[[j]] <- Zm[[j]]*NA
                for(k in 1:Km[[j]][3]) {
                    UTA[[j]][,,k] <-  upper.tri(Zm[[j]][,,1])
                } 
                UTA[[j]] <- (UTA[[j]]==1)
            }
            ## Step 2. update U
            U <- updateUm(ns, U, V, R, Zm, Km, ej, s2, eU.st, iVU, UL.Normal)
        
            ## Step 3. update V
            Vm <- updateVm(ns, U, V, Zm, Km, R, s2, eV, iVV, UTA)
            V <- Reduce(rbind, Vm)

            for(j in 1:ns){
                MU[[j]] <- M.U(list(U[[j]],U[[j]],Vm[[j]]))
                MU.state[[j]] <- M.U(list(U[[j]],U[[j]],V))
                ZU[[j]] <- ZY[[j]] - MU[[j]]
            }
            ## Step 4. update s2
            s2 <- updates2m(ns, Zm, MU, c0, d0, Km)
       
            if(constant){
                ## Step 5. update constant bhat
                bhat <- updatebm(ns, K, s, s2, B0, p, ZU)
                Zb <- Z - bhat
            }
            
            ## update hierarchical parameters
            for(j in 1:ns){
                dens.temp <- rep(NA, R)
                SS <- t(U[[j]])%*%U[[j]]
                for(r in 1:R){
                    iVU[[j]][r,r] <- 1/rgamma(1, (u0 + K[1])/2, (u1 + SS[r,r])/2)
                    dens.temp[r] <- log(dinvgamma(iVU.st[[j]][r,r], u0/2, u1/2))
                }
                ## invSS <- chol2inv(chol(SS))
                ## iVU[[j]] <- rwish(invSS, psi.U0 + Km[[j]][1])
                density.iVU.temp.holder[j] <- sum(dens.temp) ## dwishartc(U=iVU.st[[j]], nu=psi.U1,
                ##                                         S=invSS, log=TRUE)
            }
            density.iVU.holder[g] <- sum(density.iVU.temp.holder)
        
            ## hierarchical parameters for V
            for(j in 1:ns){
                Vs <- matrix(Vm[[j]], nrow=sum(ej[[j]]), ncol=R)
                SS <-  t(Vs)%*%Vs
                for(r in 1:R){
                    iVV[[j]][r,r] <- 1/rgamma(1, (v0 + Km[[j]][3])/2, (v1 + SS[r,r])/2)
                }
                ## Vs <- matrix(Vm[[j]], nrow=sum(ej[[j]]), ncol=R)
                ## SS <-  Psi.V0 +  t(Vs)%*%Vs                
                ## invSS <- chol2inv(chol(SS))
                ## iVV[[j]] <- rwish(invSS, psi.V0 + Km[[j]][3])
                eV[[j]] <- c(rMVNorm(1,apply(Vs, 2, sum)/(Km[[j]][3]+1),
                                                     solve(iVV[[j]])/(Km[[j]][3]+1)))      
            }

            ## Step 5. update s
            state.out <- updateS(iter, s, V, m, Zb, Zt, Time, 
                                  MU.state, P, s2, N.upper.tri, random.perturb)
            s <- state.out$s
            ps <- state.out$ps
            
            ## double check 
            if(length(table(s)) < ns){
                print(table(s))
                cat("Sampled s does not have all states. \n")
                s <- sort(sample(1:ns, size=K[3], replace=TRUE, prob=(rep(1, ns))))
            }
            ## Step 6. update P
            P <- updateP(s, ns, P, A0)
        }
        
        density.iVU <- log(mean(exp(density.iVU.holder)))
        if(abs(density.iVU) == Inf){
            ## cat("    Precision reinforced! \n")
            ## print(density.iVU.holder)
            density.iVU <- as.numeric(log(mean(exp(mpfr(density.iVU.holder, precBits=53)))))
        }
        if(abs(density.iVU) == Inf){
            stop("The density of a posterior ordinate reaches a computational limit and marginal likelihood computation is halted.\n")
        }

        cat("\n---------------------------------------------- \n ")
        cat("Marignal Likelihood Computation Step 2 \n")
        cat("    density.iVU: ", as.numeric(density.iVU), "\n")
        cat("---------------------------------------------- \n ")

        ##
        ## Marginal Step 3: p(eV.st|Z, eu.st, iVU.st)
        ##
        density.eV.holder <- matrix(NA, reduce.mcmc)
        for (g in 1:reduce.mcmc){
            density.eV.temp.holder <- rep(NA, ns)
            ## Step 1. update ej, Km, Zm
            for (j in 1:ns){
                ej[[j]] <- as.numeric(s==j)
                Km[[j]] <- dim(Zb[,,ej[[j]]==1])
                if(is.na(Km[[j]][3])){
                    ZY[[j]] <- array(Z[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1)))
                    Zm[[j]] <- array(Zb[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1)))
                } else{
                    ZY[[j]] <- Z[,,ej[[j]]==1]
                    Zm[[j]] <- Zb[,,ej[[j]]==1]
                }
                Km[[j]] <- dim(Zm[[j]])
                
                ## UTA array: TRUE for upper triangle
                UTA[[j]] <- Zm[[j]]*NA
                for(k in 1:Km[[j]][3]) {
                    UTA[[j]][,,k] <-  upper.tri(Zm[[j]][,,1])
                } 
                UTA[[j]] <- (UTA[[j]]==1)
            }
            ## Step 2. update U
            U <- updateUm(ns, U, V, R, Zm, Km, ej, s2, eU.st, iVU.st, UL.Normal)
        
            ## Step 3. update V
            Vm <- updateVm(ns, U, V, Zm, Km, R, s2, eV, iVV, UTA)
            V <- Reduce(rbind, Vm)

            for(j in 1:ns){
                MU[[j]] <- M.U(list(U[[j]],U[[j]],Vm[[j]]))
                MU.state[[j]] <- M.U(list(U[[j]],U[[j]],V))
                ZU[[j]] <- ZY[[j]] - MU[[j]]
            }
            ## Step 4. update s2
            s2 <- updates2m(ns, Zm, MU, c0, d0, Km)
       
            if(constant){
                ## Step 5. update constant bhat
                bhat <- updatebm(ns, K, s, s2, B0, p, ZU)
                Zb <- Z - bhat
            }
            
            ## hierarchical parameters for V
            for(j in 1:ns){
                Vs <- matrix(Vm[[j]], nrow=sum(ej[[j]]), ncol=R)
                SS <-  t(Vs)%*%Vs
                for(r in 1:R){
                    iVV[[j]][r,r] <- 1/rgamma(1, (v0 + Km[[j]][3])/2, (v1 + SS[r,r])/2)
                }
                mu.eV <- apply(Vs, 2, sum)/(Km[[j]][3]+1)
                sigma.eV <- solve(iVV[[j]])/(Km[[j]][3]+1)
                eV[[j]] <- c(rMVNorm(1, mu.eV, sigma.eV))
                density.eV.temp.holder[j] <- dmvnorm(eV.st[[j]], matrix(mu.eV, 1, R), sigma.eV, log=TRUE)
            }
            density.eV.holder[g] <- sum(density.eV.temp.holder)
            
        
            ## Step 6. update s
            state.out <- updateS(iter, s, V, m, Zb, Zt, Time, 
                                  MU.state, P, s2, N.upper.tri, random.perturb)
            ## state.out <- updateS(iter, s, V, m, Zb, Zt, Time, fast,
            ##                      MU.state, P, s2, local.type, logistic.tune, N.upper.tri, sticky)
            s <- state.out$s
            ps <- state.out$ps
            
            ## double check 
            if(length(table(s)) < ns){
                ## print(table(s))
                ## cat("Sampled s does not have all states. \n")
                s <- sort(sample(1:ns, size=K[3], replace=TRUE, prob=(rep(1, ns))))
            }
            ## Step 6. update P
            P <- updateP(s, ns, P, A0)
        }
        
        density.eV <- log(mean(exp(density.eV.holder)))
        if(abs(density.eV) == Inf){
            ## cat("    Precision reinforced! \n")
            ## print(density.eV.holder)
            density.eV <- as.numeric(log(mean(exp(mpfr(density.eV.holder, precBits=53)))))
        }
        if(abs(density.eV) == Inf){
            stop("The density of a posterior ordinate reaches a computational limit and marginal likelihood computation is halted.\n")
        }

        
        cat("\n---------------------------------------------- \n ")
        cat("Marignal Likelihood Computation Step 3 \n")
        cat("    density.eV: ", as.numeric(density.eV), "\n")
        cat("---------------------------------------------- \n ")

        ##
        ## Marginal Step 4: p(iVV.st|Z, eU.st, iVU.st, eV.st)
        ##
        density.iVV.holder <- matrix(NA, reduce.mcmc)
        for (g in 1:reduce.mcmc){
            density.iVV.temp.holder <- rep(NA, ns)
            ## Step 1. update ej, Km, Zm
            for (j in 1:ns){
                ej[[j]] <- as.numeric(s==j)
                Km[[j]] <- dim(Zb[,,ej[[j]]==1])
                if(is.na(Km[[j]][3])){
                    ZY[[j]] <- array(Z[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1)))
                    Zm[[j]] <- array(Zb[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1)))
                } else{
                    ZY[[j]] <- Z[,,ej[[j]]==1]
                    Zm[[j]] <- Zb[,,ej[[j]]==1]
                }
                Km[[j]] <- dim(Zm[[j]])
                
                ## UTA array: TRUE for upper triangle
                UTA[[j]] <- Zm[[j]]*NA
                for(k in 1:Km[[j]][3]) {
                    UTA[[j]][,,k] <-  upper.tri(Zm[[j]][,,1])
                } 
                UTA[[j]] <- (UTA[[j]]==1)
            }
            ## Step 2. update U
            U <- updateUm(ns, U, V, R, Zm, Km, ej, s2, eU.st, iVU.st, UL.Normal)
        
            ## Step 3. update V
            Vm <- updateVm(ns, U, V, Zm, Km, R, s2, eV.st, iVV, UTA)
            V <- Reduce(rbind, Vm)

            for(j in 1:ns){
                MU[[j]] <- M.U(list(U[[j]],U[[j]],Vm[[j]]))
                MU.state[[j]] <- M.U(list(U[[j]],U[[j]],V))
                ZU[[j]] <- ZY[[j]] - MU[[j]]
            }
            ## Step 4. update s2
            s2 <- updates2m(ns, Zm, MU, c0, d0, Km)
       
            if(constant){
                ## Step 5. update constant bhat
                bhat <- updatebm(ns, K, s, s2, B0, p, ZU)
                Zb <- Z - bhat
            }
            
            ## hierarchical parameters for V
            for(j in 1:ns){
                dens.temp <- rep(NA, R)
                Vs <- matrix(Vm[[j]], nrow=sum(ej[[j]]), ncol=R)
                SS <-  t(Vs)%*%Vs
                for(r in 1:R){
                    iVV[[j]][r,r] <- 1/rgamma(1, (v0 + Km[[j]][3])/2, (v1 + SS[r,r])/2)
                    dens.temp[r] <- log(dinvgamma(iVV.st[[j]][r,r], v0/2, v1/2))
                }
                ## Vs <- matrix(Vm[[j]], nrow=sum(ej[[j]]), ncol=R)
                ## SS <-  Psi.V0 +  t(Vs)%*%Vs                
                ## invSS <- chol2inv(chol(SS))
                ## iVV[[j]] <- rwish(invSS, psi.V0 + Km[[j]][3])
                ## eV[[j]] <- c(rMVNorm(1,apply(Vs, 2, sum)/(Km[[j]][3]+1),
                ## solve(iVV[[j]])/(Km[[j]][3]+1)))
                density.iVV.temp.holder[j] <- sum(dens.temp)
                ## dwishartc(U=iVV.st[[j]], nu=psi.V1, S=invSS, log=TRUE)
            }
            density.iVV.holder[g] <- sum(density.iVV.temp.holder)
           
        
            ## Step 6. update s
            state.out <- updateS(iter, s, V, m, Zb, Zt, Time, MU.state, P, s2,
                          N.upper.tri, random.perturb)
            ## state.out <- updateS(iter, s, V, m, Zb, Zt, Time, fast,
            ##                      MU.state, P, s2, local.type, logistic.tune, N.upper.tri, sticky)
            s <- state.out$s
            ps <- state.out$ps
            
            ## double check 
            if(length(table(s)) < ns){
                ## print(table(s))
                ## cat("Sampled s does not have all states. \n")
                s <- sort(sample(1:ns, size=K[3], replace=TRUE, prob=(rep(1, ns))))
            }
            ## Step 6. update P
            P <- updateP(s, ns, P, A0)
        }
        
        density.iVV <- log(mean(exp(density.iVV.holder)))
        if(abs(density.iVV) == Inf){
            ## cat("    Precision reinforced! \n")
            ## print(density.iVV.holder)
            density.iVV <- as.numeric(log(mean(exp(mpfr(density.iVV.holder, precBits=53)))))
        }
        if(abs(density.iVV) == Inf){
            stop("The density of a posterior ordinate reaches a computational limit and marginal likelihood computation is halted.\n")
        }


        cat("\n---------------------------------------------- \n ")
        cat("Marignal Likelihood Computation Step 4 \n")
        cat("    density.iVV: ", as.numeric(density.iVV), "\n")
        cat("---------------------------------------------- \n ")

        if(constant){
            ##
            ## Marginal Step 5: p(bhat.st|Z, eu.st, iVU.st, eV.st, iVV.st)
            ##
            density.bhat.holder <- matrix(NA, reduce.mcmc)
            for (g in 1:reduce.mcmc){
                ## Step 1. update ej, Km, Zm
                for (j in 1:ns){
                    ej[[j]] <- as.numeric(s==j)
                    Km[[j]] <- dim(Zb[,,ej[[j]]==1])
                    if(is.na(Km[[j]][3])){
                        ZY[[j]] <- array(Z[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1)))
                        Zm[[j]] <- array(Zb[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1)))
                    } else{
                        ZY[[j]] <- Z[,,ej[[j]]==1]
                        Zm[[j]] <- Zb[,,ej[[j]]==1]
                    }
                    Km[[j]] <- dim(Zm[[j]])
                    
                    ## UTA array: TRUE for upper triangle
                    UTA[[j]] <- Zm[[j]]*NA
                    for(k in 1:Km[[j]][3]) {
                        UTA[[j]][,,k] <-  upper.tri(Zm[[j]][,,1])
                    } 
                    UTA[[j]] <- (UTA[[j]]==1)
                }
                ## Step 2. update U
                U <- updateUm(ns, U, V, R, Zm, Km, ej, s2, eU.st, iVU.st, UL.Normal)
                
                ## Step 3. update V
                Vm <- updateVm(ns, U, V, Zm, Km, R, s2, eV.st, iVV.st, UTA)
                V <- Reduce(rbind, Vm)
                
                for(j in 1:ns){
                    MU[[j]] <- M.U(list(U[[j]],U[[j]],Vm[[j]]))
                    MU.state[[j]] <- M.U(list(U[[j]],U[[j]],V))
                ZU[[j]] <- ZY[[j]] - MU[[j]]
                }
                ## Step 4. update s2
                s2 <- updates2m(ns, Zm, MU, c0, d0, Km)
                
                ## Step 5. update constant bhat
                cV <- 1/(sum(sapply(1:ns, function(j){prod(K[1:2])*length(s == j)/s2[j]})) +  1/B0)
                cE <- cV*sum(unlist(lapply(ZU, sum)))
                bhat <- rnorm(1, cE, sqrt(cV)) 
                Zb <- Z - bhat
                density.bhat.holder[g] <- dnorm(bhat.st, cE, sqrt(cV), log=TRUE)
                
                ## Step 6. update s
                state.out <- updateS(iter, s, V, m, Zb, Zt, Time, MU.state, P, s2,
                                     N.upper.tri, random.perturb)
                ## state.out <- updateS(iter, s, V, m, Zb, Zt, Time, fast,
                ##              MU.state, P, s2, local.type, logistic.tune, N.upper.tri, sticky)
                s <- state.out$s
                ps <- state.out$ps
                
                ## double check 
                if(length(table(s)) < ns){
                    ## print(table(s))
                    ## cat("Sampled s does not have all states. \n")
                    s <- sort(sample(1:ns, size=K[3], replace=TRUE, prob=(rep(1, ns))))
                }
                ## Step 6. update P
                P <- updateP(s, ns, P, A0)
            }
            
            density.bhat <- log(mean(exp(density.bhat.holder)))
            if(abs(density.bhat) == Inf){
                ## cat("    Precision reinforced! \n")
                ## print(density.bhat.holder)
                density.bhat <- as.numeric(log(mean(exp(mpfr(density.bhat.holder, precBits=53)))))
            }
            if(abs(density.bhat) == Inf){
                stop("The density of a posterior ordinate reaches a computational limit and marginal likelihood computation is halted.\n")
            }
            
            
            cat("\n---------------------------------------------- \n ")
            cat("Marignal Likelihood Computation Step 5 \n")
            cat("    density.bhat: ", as.numeric(density.bhat), "\n")
            cat("---------------------------------------------- \n ")
        }
        
        ##
        ## Marginal Step 6: p(s2.st|Z, eu.st, iVU.st, eV.st, iVV.st, bhat.st)
        ##
        density.Sigma.holder <- matrix(NA, reduce.mcmc)
        for (g in 1:reduce.mcmc){
            ## Step 1. update ej, Km, Zm
            for (j in 1:ns){
                ej[[j]] <- as.numeric(s==j)
                Km[[j]] <- dim(Zb[,,ej[[j]]==1])
                if(is.na(Km[[j]][3])){
                    ZY[[j]] <- array(Z[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1)))
                    Zm[[j]] <- array(Zb[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1)))
                } else{
                    ZY[[j]] <- Z[,,ej[[j]]==1]
                    Zm[[j]] <- Zb[,,ej[[j]]==1]
                }
                Km[[j]] <- dim(Zm[[j]])
                
                ## UTA array: TRUE for upper triangle
                UTA[[j]] <- Zm[[j]]*NA
                for(k in 1:Km[[j]][3]) {
                    UTA[[j]][,,k] <-  upper.tri(Zm[[j]][,,1])
                } 
                UTA[[j]] <- (UTA[[j]]==1)
            }
            ## Step 2. update U
            U <- updateUm(ns, U, V, R, Zm, Km, ej, s2, eU.st, iVU.st, UL.Normal)
        
            ## Step 3. update V
            Vm <- updateVm(ns, U, V, Zm, Km, R, s2, eV.st, iVV.st, UTA)
            V <- Reduce(rbind, Vm)

            for(j in 1:ns){
                MU[[j]] <- M.U(list(U[[j]],U[[j]],Vm[[j]]))
                MU.state[[j]] <- M.U(list(U[[j]],U[[j]],V))
                ZU[[j]] <- ZY[[j]] - MU[[j]]
            }
            ## Step 4. update s2
            ZEE <- as.list(rep(NA, ns))
            density.Sigma <- rep(NA, ns)
            for(j in 1:ns){
                ZEE[[j]] <- Zm[[j]] - MU[[j]]
                EE <- c(ZEE[[j]])
                shape <- (c0+prod(Km[[j]]))/2
                scale <- (d0+ sum(EE^2))/2
                s2[j] <- 1/rgamma(1, shape, scale)
                density.Sigma[j] <- log(dinvgamma(s2.st[j], shape, scale))
            }  
            density.Sigma.holder[g] <- sum(density.Sigma)
       
            ## Step 5. update constant bhat
            ## Zb.st <- Z - bhat.st
        
            ## Step 6. update s
            state.out <- updateS(iter, s, V, m, Zb.st, Zt, Time, MU.state, P, s2,
                          N.upper.tri, random.perturb)
            ## state.out <- updateS(iter, s, V, m, Zb.st, Zt, Time, fast,
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
            P <- updateP(s, ns, P, A0)
        }
        
        density.Sigma <- log(mean(exp(density.Sigma.holder)))
        if(abs(density.Sigma) == Inf){
            ## cat("    Precision reinforced! \n")
            ## print(density.Sigma.holder)
            density.Sigma <- as.numeric(log(mean(exp(mpfr(density.Sigma.holder, precBits=53)))))
        }
        if(abs(density.Sigma) == Inf){
            stop("The density of a posterior ordinate reaches a computational limit and marginal likelihood computation is halted.\n")
        }
        
        cat("\n---------------------------------------------- \n ")
        cat("Marignal Likelihood Computation Step 6 \n")
        cat("    density.Sigma: ", as.numeric(density.Sigma), "\n")
        cat("---------------------------------------------- \n ")

        
        ##
        ## Marginal Step 7: p(P.st| Z, eU.st, iVU.st, eV.st, iVV.st, b.st, s2.st)
        ##
        density.P.holder <- matrix(NA, reduce.mcmc)
        for (g in 1:reduce.mcmc){
            density.P.temp.holder <- rep(NA, ns)
            ## Step 1. update ej, Km, Zm
            for (j in 1:ns){
                ej[[j]] <- as.numeric(s==j)
                Km[[j]] <- dim(Zb[,,ej[[j]]==1])
                if(is.na(Km[[j]][3])){
                    ZY[[j]] <- array(Z[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1)))
                    Zm[[j]] <- array(Zb[,,ej[[j]]==1], dim = c(Km[[j]][1], Km[[j]][2], sum(ej[[j]]==1)))
                } else{
                    ZY[[j]] <- Z[,,ej[[j]]==1]
                    Zm[[j]] <- Zb[,,ej[[j]]==1]
                }
                Km[[j]] <- dim(Zm[[j]])
                
                ## UTA array: TRUE for upper triangle
                UTA[[j]] <- Zm[[j]]*NA
                for(k in 1:Km[[j]][3]) {
                    UTA[[j]][,,k] <-  upper.tri(Zm[[j]][,,1])
                } 
                UTA[[j]] <- (UTA[[j]]==1)
            }
            ## Step 2. update U
            U <- updateUm(ns, U, V, R, Zm, Km, ej, s2.st, eU.st, iVU.st, UL.Normal)
        
            ## Step 3. update V
            Vm <- updateVm(ns, U, V, Zm, Km, R, s2.st, eV.st, iVV.st, UTA)
            V <- Reduce(rbind, Vm)

            for(j in 1:ns){
                MU[[j]] <- M.U(list(U[[j]],U[[j]],Vm[[j]]))
                MU.state[[j]] <- M.U(list(U[[j]],U[[j]],V))
                ZU[[j]] <- ZY[[j]] - MU[[j]]
            }
            ## Step 4. update s2
       
            ## Step 5. update constant bhat
        
            ## Step 6. update s
            state.out <- updateS(iter, s, V, m, Zb.st, Zt, Time, MU.state, P, s2.st,
                          N.upper.tri, random.perturb)
            ##            state.out <- updateS(iter, s, V, m, Zb.st, Zt, Time, fast,
            ##            MU.state, P, s2.st, local.type, logistic.tune, N.upper.tri, sticky)
            s <- state.out$s
            ps <- state.out$ps
            
            ## double check 
            if(length(table(s)) < ns){
                ## print(table(s))
                ## cat("Sampled s does not have all states. \n")
                s <- sort(sample(1:ns, size=K[3], replace=TRUE, prob=(rep(1, ns))))
            }
            ## Step 6. update P
            P <- updateP(s, ns, P, A0)
            swit  <-  switchg(s) 
            for (j in 1:ns){
                swit1 <-  A0[j,] + swit[j,]
                pj <-  rdirichlet.cp(1, swit1)
                P[j,] <-  pj
                if(j < ns){
                    shape1 <-  swit1[j]
                    shape2 <-  swit1[j + 1]   
                    density.P.temp.holder[j] <- log(dbeta(P.st[j], shape1, shape2))
                }
            }
            density.P.holder[g] <- sum(density.P.temp.holder, na.rm=TRUE)
            
        }
        
        density.P <- log(mean(exp(density.P.holder)))
        if(abs(density.P) == Inf){
            ## cat("    Precision reinforced! \n")
            ## print(density.P.holder)
            density.P <- as.numeric(log(mean(exp(mpfr(density.P.holder, precBits=53)))))
        }
         if(abs(density.P) == Inf){
            stop("The density of a posterior ordinate reaches a computational limit and marginal likelihood computation is halted.\n")
        }

        cat("\n---------------------------------------------- \n ")
        cat("Marignal Likelihood Computation Step 7 \n")
        cat("    density.P: ", as.numeric(density.P), "\n")
        cat("---------------------------------------------- \n ")

        ## Prior ordinate
        iVV0 <- iVU0 <- diag(R) ; eV0 <- eU0 <- rep(0,R) ## ; psi0 <- rep(1,R) ; Psi0 <- R+1
        if (is.null(a) | is.null(b)) {
            expected.duration <- round(K[3]/(m + 1))
            b <- 0.1
            a <- b * expected.duration
        }
        density.eU.prior <- density.eV.prior <- density.iVU.prior <- density.iVV.prior <- Sigma.prior <- rep(NA, ns)
        density.P.prior <- rep(NA, ns-1)
        for (j in 1:ns){
            Km.st[[j]][3] <- ifelse(is.na(Km[[j]][3]), 1, Km.st[[j]][3])
            ## Prior ordinate
            density.eU.prior[j] <- log(dmvnorm(eU.st[[j]], eU0, iVU0)) 
            density.eV.prior[j] <- log(dmvnorm(eV.st[[j]], eV0, iVV0))

            dens.temp.iVU <- dens.temp.iVV <- rep(NA, R)
            for(r in 1:R){
                dens.temp.iVU[r] <- log(dinvgamma(iVU.st[[j]][r,r], u0/2, u1/2))
                dens.temp.iVV[r] <- log(dinvgamma(iVV.st[[j]][r,r], v0/2, v1/2))
            }
            density.iVU.prior[j] <- sum(dens.temp.iVU)
            density.iVV.prior[j] <- sum(dens.temp.iVV)
            Sigma.prior[j] <- log(dinvgamma(s.st[j], c0/2, (d0)/2))
            if(j < ns){
                density.P.prior[j] <- log(dbeta(P.st[j], a, b)) ## p = 1
            }
        }
        if(constant){
            density.beta.prior <- dnorm(bhat.st, b0, sqrt(B0), log=TRUE) ## p = 1
        }
        density.Sigma.prior <- sum(Sigma.prior)

        ## marginal prior density
        if(constant){
            logprior <- sum(density.eU.prior) + sum(density.eV.prior) + sum(density.iVU.prior) +
                sum(density.iVV.prior) +
                density.beta.prior + density.Sigma.prior + sum(density.P.prior)
            logdenom <- (density.eU + density.eV + density.iVU + density.iVV +
                         density.bhat + density.Sigma + density.P)
        } else{
            logprior <- sum(density.eU.prior) + sum(density.eV.prior) + sum(density.iVU.prior) +
                sum(density.iVV.prior) +
                density.Sigma.prior + sum(density.P.prior)
            logdenom <- (density.eU + density.eV + density.iVU + density.iVV +
                         density.Sigma + density.P)
        }
        ## logmarglike <- (loglike + logprior) - logdenom;
        logmarglike.upper <- (loglike.upper + logprior) - logdenom
        
        cat("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n")
        ## cat("    log marginal likelihood = (loglike + logprior) - (density.parameters) \n")
        ## cat("    log marginal likelihood : ", as.numeric(logmarglike), "\n")
        cat("    log marginal likelihood : ", as.numeric(logmarglike.upper), "\n")
        cat("    logprior : ", as.numeric(logprior), "\n")
        cat("    log posterior density : ", as.numeric(logdenom), "\n")
        cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n")
                
    }
    #########################################################################
    ## flip V and U based on the size of mse(V)
    ## so that the first principal axis comes at the first column of V and U
    ## In that way, the interpretation of cp results makes sense.
    #########################################################################
    output <- abind(MU.st) ## lapply(MU.record, function(x){x/nss}) ## MU.record/nss
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
    attr(output, "DIC") <- sum(as.numeric(Z.DIC))
    attr(output, "Waic.out") <- Waic.out
    attr(output, "prob.state") <- ps.store/nss
    ## attr(output, "loglike") <- loglike
    attr(output, "loglike") <- loglike.upper
    ## attr(output, "logmarglike") <- logmarglike
    attr(output, "logmarglike") <- logmarglike.upper 
    ## cat("elapsed time for m = ", m, " is ", proc.time() - ptm, "\n")
    return(output)
}

