#' Degree-corrected multilinear tensor model
#'
#' NetworkStatic implements a degree-corrected Bayesian multilinear tensor decomposition method

#' @param Y Reponse tensor 
#' @param R Dimension of latent space. The default is 2. 
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

#' @param marginal If \code{marginal = TRUE}, the log marignal likelihood is computed using the method of Chib (1995).
    
#' @param DIC If \code{DIC = TRUE}, the deviation information criterion is computed.
    
#' @param Waic If \code{Waic = TRUE}, the Watanabe information criterion is computed.
#'
#' @return An mcmc object that contains the posterior sample. This object can
#' be summarized by functions provided by the coda package.  The object
#' contains an attribute \code{Waic.out} that contains results of WAIC and the log-marginal
#' likelihood of the model (\code{logmarglike}).
#'
#' @seealso \code{\link{NetworkChange}}
#'
#' @references  Jong Hee Park and Yunkyun Sohn. 2017. "Detecting Structural Change
#' in Network Time Series Data using Bayesian Inference." Working Paper.
#'
#' Peter D. Hoff 2011. "Hierarchical Multilinear Models for Multiway Data."
#' \emph{Computational Statistics \& Data Analysis}. 55: 530-543.
#'
#' Sumio Watanabe. 2010. "Asymptotic equivalence of Bayes cross validation and widely
#' applicable information criterion in singular learning theory."
#' \emph{Journal of Machine Learning Research}. 11: 3571-3594.

#' Siddhartha Chib. 1995. ``Marginal Likelihood from the Gibbs Output.''
#' \emph{Journal of the American Statistical Association}. 90: 1313-1321.

#' @export
#' @examples
#'
#'    \dontrun{
#'    set.seed(1973)
#'
#'    ## generate an array with three constant blocks
#'    Y <- MakeBlockNetworkChange(n=10, shape=10, T=10, type ="constant")
#'    G <- 100 ## Small mcmc scans to save time
#'    out0 <- NetworkStatic(Y, R=2, mcmc=G, burnin=G, verbose=G)
#'
#'    ## recovered latent blocks
#'    Kmeans(out0, n.cluster=3, main="Recovered Blocks")
#' 
#'    ## contour plot of latent node positions
#'    plotContour(out0)
#'
#'    ## plot latent node positions
#'    plotU(out0)
#'
#'    ## plot layer-specific network connection rules
#'    plotV(out0)
#'    }
#'
NetworkStatic <- function(Y, R=2, mcmc = 100, burnin = 100, verbose = 0,thin = 1, 
                          degree.normal="eigen", UL.Normal = "Orthonormal",
                          plotUU = FALSE, plotZ = FALSE, constant=FALSE, 
                          b0 = 0, B0 = 1, c0 = NULL, d0 = NULL,
                          u0 = NULL, u1 = NULL, v0 = NULL, v1 = NULL,
                          marginal = FALSE, DIC = FALSE,  Waic=FALSE){
    
    ## time keeper
    ## ptm <- proc.time()

    ##
    ## function call
    ##
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    
    ##
    ## MCMC controllers
    ##
    totiter <- mcmc + burnin
    nstore <- mcmc/thin    
    reduce.mcmc <- nstore
    if(is.na(dim(Y)[3])){
        Y <- array(Y, dim=c(dim(Y)[1], dim(Y)[2], 1))
    }
    K <- dim(Y)
    n <- dim(Y)[1] 
    Z <- Y
    MU.record <- Y*0
    nss <- 0


    if(R==1 & UL.Normal == "Orthonormal"|| R==1 & UL.Normal == "Normal"){
        stop("If R=1, please set UL.Normal=FALSE.")
    }
    ##
    ## Degree normalization
    ##
    if(degree.normal == "eigen"){
        for(k in 1:K[3]){
            ee <- eigen(Y[,,k])
            Z[,,k] <- Y[,,k] - ee$values[1] * outer(ee$vectors[,1], ee$vectors[,1])
        }
    }
    ## if Lsym
    if(degree.normal == "Lsym"){
        for(k in 1:K[3]){
            Yk <-as.matrix(Y[,,k])
            V <- colSums(Yk)
            L <- diag(V) - Yk
            Z[,,k] <- diag(V^(-.5))%*% L %*%diag(V^(-.5))
        }
    }
    ## if Modul
    gamma.par = 1
    if(degree.normal == "Modul"){
        for(k in 1:K[3]){
          Yk <- as.matrix(Y[,,k])
          yk <- as.vector(apply(Yk, 2, sum))
          ym <- sum(yk)
          Z[,,k] <- Yk - gamma.par*(yk%o%yk)/ym
      }
    }
    
    ##
    ## eigen decomposition for initial values of U, V, MU
    ##
    out <- startUV(Z, R, K)
    U <- out[[1]]
    V <- out[[2]]
    MU <- M.U(list(U,U,V))

    ## unique values of Y
    uy <- sort(unique(c(Y)))
    
    ## UTA array: TRUE for upper triangle
    UTA <- Z*NA
    for(k in 1:K[3]) {
        UTA[,,k] <-  upper.tri(Z[,,1] )
    } 
    UTA <- (UTA==1)
      
    ##
    ## MCMC holders
    ##
    Umat <- matrix(NA, nstore, dim(Y)[[1]]*R)
    Vmat <- matrix(NA, nstore, R*dim(Y)[[3]])
    s2mat <- bmat <- matrix(NA, nstore)
    eUmat <- matrix(NA, nstore, R)
    eVmat <- matrix(NA, nstore, R)
    iVUmat <- matrix(NA, nstore, R*R)
    iVVmat <- matrix(NA, nstore, R*R)

    ##
    ## prior
    ##
    if(is.null(c0)){
        c0 <- 1
    }
    if(is.null(d0)){
        d0 <- var(as.vector(Z - MU))
    }
    if(is.null(u0)){
        u0 <- 10
    }
    if(is.null(u1)){
        u1 <- 1 
    }
    if(is.null(v0)){
        v0 <- 4
    }
    if(is.null(v1)){
        v1 <- 2
    }
    
    ## initialize parameters
    if(constant){
        bhat <- mean(c(Z))
        Zb <- Z - bhat
    }else{
        bhat = 0
        Zb <- Z
    }
    X <- array(1, dim=c(K, 1))
    p <- dim(X)[4]
    XtX <- prod(K) 
    rm(X)
    s2 <- d0
    iVV <- iVU <- diag(R) ; eV <- eU <- rep(0,R) ; 
    
    ##
    ## Model diagnositics
    ##
    Z.loglike <- NA
    if(Waic == TRUE){
        N.upper.tri <- K[1]*(K[1]-1)/2
        Z.loglike <- matrix(NA, nstore, N.upper.tri*K[3])
    }
    UTAsingle <-  upper.tri(Z[,,1])          

    if(verbose !=0){
        cat("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n")
        cat("\t NetworkStatic MCMC Sampler Starts! \n")
        ## cat("\t function called: ")
        ## print(call)
        cat("\t degree normalization: ", degree.normal, "\n")
        cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n")
    }
    ## ----------------------------------------------
    ## MCMC loop starts!  
    ## ----------------------------------------------
    for(iter in 1:totiter) {   

        ## Step 1: update bhat
        if(constant){
            bhat <- updateb(Z, MU, s2, XtX, b0, B0)
            Zb <- Z - bhat
        }
        
        ## Step 2: update U
        U <- updateU(K, U, V, R, Zb, s2, eU, iVU)
        if (UL.Normal == "Normal"){
            U <- Unormal(U)
        }else if(UL.Normal == "Orthonormal"){
            U <- GramSchmidt(U)
        }
        
        ## Step 3: update V
        V <- updateV(Zb, U, R, K, s2, eV, iVV, UTA)         
        ## update MU
        MU <- M.U(list(U,U,V))
        
        ## Step 4: update eU and iVU
        SS <-  t(U)%*%U
        for(r in 1:R){
            iVU[r,r] <- 1/rgamma(1, (u0 + K[1])/2, (u1+ SS[r,r])/2)
        }
        ## invSS <- chol2inv(chol(SS))
        ## iVU <- rwish(invSS, psi.U1)
        eU <- c(rMVNorm(1, apply(U,2,sum)/(K[1]+1), solve(iVU)/(K[1]+1)))
        
        ## Step 5: update eV and iVV
        SS <-  t(V)%*%V ## (K[3]-1)*cov(V)   +K[3]*msi/(K[3]+1)
        for(r in 1:R){
            iVV[r,r] <- 1/rgamma(1, (v0 + K[3])/2, (v1 + SS[r,r])/2)
        }
        ## invSS <- chol2inv(chol(SS))
        ## iVV <- rwish(invSS, psi.V1)
        eV <- c(rMVNorm(1,apply(V,2,sum)/(K[3]+1), solve(iVV)/(K[3]+1)))
        
        ## Step 6: update s2
        ZE <- Zb - MU 
        s2 <- rinvgamma(1, (c0 + XtX)/2, (d0+ sum(c(ZE)^2))/2  )
        
        ## report
        if (verbose!= 0 &iter %% verbose == 0){
            cat("\n----------------------------------------------",'\n')
            cat("    iteration = ", iter, '\n')
            if(constant){
                cat("    beta = ", bhat,'\n')
            }
            cat("    sigma2 = ", s2 ,'\n')
            if(plotZ == TRUE & plotUU == TRUE){
                par(mfrow=c(1, 2))
            }
            if(plotZ == TRUE){
                plot(density(c(Z)), lwd=2, main="Density of Z and MU")
                lines(density(c(MU)), col="red")
                legend("topright", legend=c("Z", "MU"), col=c("black", "red"), lty=1, lwd=1)
            }
            if(plotUU == TRUE){
                plot(U[,1], U[, 2], pch=19, cex=1); abline(v=0, col=2); abline(h=0, col=2)
            }
            ## cat("V(1,1,...R,R) = ", V, '\n')
            ## dens <- sum(log(dnorm(c(Z),c(MU),sqrt(s2))) )
            ## cat("log density of (Z|mu, sigma) = ", dens, '\n')
            cat("----------------------------------------------",'\n')
        }
        
        ## save
        if (iter > burnin & (iter-burnin)%%thin == 0){
            nss <- nss+1
            MU.record <- MU.record + MU
            if(constant){
                bmat[(iter-burnin)/thin] <- bhat
            }
            Umat[(iter-burnin)/thin, ] <- as.vector(U)
            Vmat[(iter-burnin)/thin, ] <- as.vector(V)
            eUmat[(iter-burnin)/thin, ] <- eU
            eVmat[(iter-burnin)/thin, ] <- eV
            iVUmat[(iter-burnin)/thin, ] <- iVU
            iVVmat[(iter-burnin)/thin, ] <- iVV
            
            
            s2mat[(iter-burnin)/thin] <- s2
            if(Waic == TRUE){
                d <- sapply(1:K[3], function(t){dnorm(c(Zb[,,t][UTAsingle]),
                                                      mean = c(MU[,,t][UTAsingle]),
                                                      sd=sqrt(s2), log=TRUE)})
                ## d <- sapply(1:K[3], function(t){log(dnorm(c(Zb1[,,t]), mean = c(MUU[,,t]), sqrt(s2)))})
                Z.loglike[(iter-burnin)/thin, ] <- d
            }
        }
    }## end of MCMC loop

##########################################################################
    ## Likelihood computation
##########################################################################
    logmarglike.upper <- loglike.upper <- NA
    
    ## Prepare Stars
    if(constant){
        bhat.st <- mean(bmat)
        Zb.st <- Z - bhat.st
    }else{
        Zb.st <- Z
    }
    Sigma.st <- mean(s2mat)
    eU.st <-apply(eUmat, 2, mean)
    eV.st <- apply(eVmat, 2, mean)
    iVU.st <- matrix(apply(iVUmat, 2, mean), R, R)
    iVV.st <- matrix(apply(iVVmat, 2, mean), R, R)
    U.st <- matrix(apply(Umat, 2, mean), K[1], R)
    L.st <- matrix(apply(Vmat, 2, mean), K[3], R)
    ## if (UL.Normal == "Normal"){
    ##     U.st <- Unormal(U.st)
    ## }else if(UL.Normal == "Orthonormal"){
    ##     U.st <- GramSchmidt(U.st)
    ## }
    MU.st <- M.U(list(U.st, U.st, L.st))
    EE.st <- c(Zb.st - MU.st) 
    
    
    ## compute marginal likelihood
    loglike.upper <- sum(sapply(1:K[3], function(t){dnorm(c(Zb.st[,,t][UTAsingle]),
                                                          mean = c(MU.st[,,t][UTAsingle]),
                                                          sd=sqrt(Sigma.st), log=TRUE)}))
    ## loglike <- sum(sapply(1:K[3], function(t){dnorm(c(Zb.st[,,t]),
    ##                                                mean = c(MU.st[,,t]), sqrt(Sigma.st), log=TRUE)}))
    ## cat("    loglike: ", as.numeric(loglike), "\n")
    if(verbose !=0){
        cat("    loglike : ", as.numeric(loglike.upper), "\n")
    }
    
##########################################################################
    ## DIC and Waic
##########################################################################
    ## drop some objects not used in the following   
    Waic.out <- Z.DIC <- NA
    if (DIC == TRUE){
        ## DIC computation = 2*D_average - D_theta^hat
        ## D_theta^hat: posterior estimates
        ## M.mean <- MU.record/nss
        ## s2.mean <- mean(s2mat) ## mean((trueY - M.mean)^2 )
        Z.D.hat <- loglike.upper ## sum(sapply(1:K[3], function(t){log(mean(dnorm(c(Zb[,,t]),
                                 ##                                mean = c(M.mean[,,t]), sqrt(s2.mean))))}))
        Z.D.average <- -2*sum(colMeans(Z.loglike))
        Z.DIC <- 2*Z.D.average - Z.D.hat
        
        cat("----------------------------------------------",'\n')
        ##if(DIC.report==TRUE){
        ## cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n")
        ## cat("log p(y|theta.Bayes): ", sum(as.numeric(logpy)), "\n")
        ## cat("Deviance Information = -2*log p(trueY|theta) \n")
        ## cat("Deviance Information = -2*log p(Z|theta) \n")
        ## cat("D.average: ", as.numeric(Z.D.average), "\n")
        ## cat("D.hat: ", as.numeric(Z.D.hat), "\n")
        cat("    DIC: ", as.numeric(Z.DIC), "\n")
        ## cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n")
        ## }
    }
    if(Waic == TRUE){
        ## Waic computation
        Waic.out <- waic(Z.loglike)$total
        rm(Z.loglike)
        
        ## 
        ## cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n")
        if(verbose !=0){
            cat("    Waic: ", Waic.out[1], "\n")
            ## cat("lpd: ", Waic.out[3], "\n")
            ## cat("p_Waic: ", Waic.out[4], "\n")
            ## cat("elapsed time for m = 0 is ", proc.time() - ptm, "\n")
            ## cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n")
            cat("----------------------------------------------",'\n')
        }
    }
    
    output <- MU.record/nss
    names(output) <- "MU"
    ## rm(MU.record) 

##########################################################################
    ## Marginal Likelihood computation
##########################################################################

    if(marginal == TRUE){
       
        
        ## Hierarchical parameters must be properly considered.
        ## 1. p(eu.st|Z)
        ## 2. p(iVU.st|Z, eu.st) = p()
        ## 3. p(ed.st|Z, eu.st, iVU.st)
        ## 4. p(iVD.st|Z, eu.st, iVU.st, ed.st)
        ## 5. same
        ##
        ## Marginal Step 1: p(eU.st|Z)
        ##
        density.eU.holder <- matrix(NA, reduce.mcmc)
        for (g in 1:reduce.mcmc){
            ## if(constant){
            ##     Zbg <- Z - bmat[g]
            ## }
            Vg <- matrix(Vmat[g, ], K[3], R)
            Ug <- matrix(Umat[g,], K[1], R)

            ## Evaluate eU.st
            SS <-  t(Ug) %*% Ug ## (K[1]-1)*cov(Ug) + K[1]*msi/(K[1]+1)
            ## invSS <- chol2inv(chol(SS))
            ## iVU <- rwish(invSS, psi.U1)
            for(r in 1:R){
                iVU[r,r] <- 1/rgamma(1, (u0 + K[1])/2, (u1+ SS[r,r])/2)
            }
            mu <- matrix(apply(Ug,2,sum)/(K[1]+1), R, 1)
            sigma <- solve(iVU)/(K[1]+1)       
            ## if(Precision == TRUE){    
            ##     density.eU.holder[g] <- dmvnorm_log(x = matrix(eU.st, R, 1),
            ##                                         mu, sigma)
            ## }else{
            density.eU.holder[g] <- dmvnorm(eU.st, matrix(mu, 1, R), sigma, log=TRUE)
            ## }
        }
        density.eU <- log(mean(exp(density.eU.holder)))
        if(abs(density.eU) == Inf){
            cat("    Precision reinforced! \n")
            print(density.eU.holder)
            density.eU <- as.numeric(log(mean(exp(mpfr(density.eU.holder, precBits=53)))))
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
            ## Step 1: update bhat
            if(constant){
                bhat <- updateb(Z, MU, s2, XtX, b0, B0)
                Zb <- Z - bhat
            }
            
            ## Step 2: update U
            U <- updateU(K, U, V, R, Zb, s2, eU.st, iVU)
            if (UL.Normal == "Normal"){
                U <- Unormal(U)
            }else if(UL.Normal == "Orthonormal"){
                U <- GramSchmidt(U)
            }
            
            ## Step 3: update V
            V <- updateV(Zb, U, R, K, s2, eV, iVV, UTA)         
            ## update MU
            MU <- M.U(list(U,U,V))
            
            ## Step 4: update eU and iVU
            SS <- t(U)%*%U
            dens.temp <- rep(NA, R)
            for(r in 1:R){
                iVU[r,r] <- 1/rgamma(1, (u0 + K[1])/2, (u1+ SS[r,r])/2)
                dens.temp[r] <- log(dinvgamma(iVU.st[r,r], u0/2, u1/2))
            }
                ## invSS <- chol2inv(chol(SS))
                ## iVU[[j]] <- rwish(invSS, psi.U0 + Km[[j]][1])
            density.iVU.holder[g] <- sum(dens.temp)

            
            ## invSS <- chol2inv(chol(SS))
            ## iVU <- rwish(invSS, psi.U1)
            ## density.iVU.holder[g] <- dwishartc(U=iVU.st, nu=psi.U1, S=invSS, log=TRUE) ## log(dwish(iVU.st, v=psi.U1, S=invSS))

            ## dwishart(Omega = iVU.st, nu = psi0 + K[1], S=SS)
            ## eU <- c(rMVNorm(1, apply(U,2,sum)/(K[1]+1), solve(iVU)/(K[1]+1)))
            
            ## Step 5: update eV and iVV
            SS <-  t(V)%*%V ## (K[3]-1)*cov(V)   +K[3]*msi/(K[3]+1)
            for(r in 1:R){
                iVV[r,r] <- 1/rgamma(1, (v0 + K[3])/2, (v1+ SS[r,r])/2)
            }
            ## invSS <- chol2inv(chol(SS))
            ## iVV <- rwish(invSS, psi.V1)
            eV <- c(rMVNorm(1,apply(V,2,sum)/(K[3]+1), solve(iVV)/(K[3]+1)))
            
            ## Step 6: update s2
            ZE <- Zb - MU 
            s2 <- rinvgamma(1, (c0 + XtX)/2, (d0+ sum(c(ZE)^2))/2  )        
        }
        
        density.iVU <- log(mean(exp(density.iVU.holder)))
        if(abs(density.iVU) == Inf){
            cat("    Precision reinforced! \n")
            print(density.iVU.holder)
            density.iVU <- as.numeric(log(mean(exp(mpfr(density.iVU.holder, precBits=53)))))
        }
        
        cat("\n---------------------------------------------- \n ")
        cat("Marignal Likelihood Computation Step 2 \n")
        cat("    density.iVU: ", as.numeric(density.iVU), "\n")
        cat("---------------------------------------------- \n ")


        ##
        ## Marginal Step 3: p(eV.st|Z, eU.st, iVU.st)
        ##
        density.eV.holder <- matrix(NA, reduce.mcmc)
        for (g in 1:reduce.mcmc){
            
            ## Step 1: update bhat
            if(constant){
                bhat <- updateb(Z, MU, s2, XtX, b0, B0)
                Zb <- Z - bhat
            }
            
            ## Step 2: update U
            U <- updateU(K, U, V, R, Zb, s2, eU.st, iVU.st)
            if (UL.Normal == "Normal"){
                U <- Unormal(U)
            }else if(UL.Normal == "Orthonormal"){
                U <- GramSchmidt(U)
            }
            
            ## Step 3: update V
            V <- updateV(Zb, U, R, K, s2, eV, iVV, UTA)         
            ## update MU
            MU <- M.U(list(U,U,V))
            
            ## Step 4: update eU and iVU
            ## SS <-  diag(psi0, nrow=R) +  t(U)%*%U 
            ## iVU <- rwish(solve(SS), Psi0+K[1])
            ## eU <- c(rMVNorm(1, apply(U,2,sum)/(K[1]+1), solve(iVU)/(K[1]+1)))
            
            ## Step 5: update eV and iVV
            SS <-  t(V)%*%V ## (K[3]-1)*cov(V)   +K[3]*msi/(K[3]+1)
            ## invSS <- chol2inv(chol(SS))
            ## iVV <- rwish(invSS, psi.V1)
            for(r in 1:R){
                iVV[r,r] <- 1/rgamma(1, (v0 + K[3])/2, (v1+ SS[r,r])/2)
            }
            mu.eV <- matrix(apply(V,2,sum)/(K[3]+1), R, 1)
            sigma.eV <- solve(iVV)/(K[3]+1)
            eV <- c(rMVNorm(1, mu.eV, sigma.eV))
            ## if(Precision == TRUE){                
            ##     density.eV.holder[g] <- dmvnorm_log(eV.st, mu.eV, sigma.eV)
            ## }else{
            density.eV.holder[g] <- dmvnorm(eV.st, matrix(mu.eV, 1, R), sigma.eV, log=TRUE)
            ##  }
             
            ## Step 6: update s2
            ZE <- Zb - MU 
            s2 <- rinvgamma(1, (c0 + XtX)/2, (d0+ sum(c(ZE)^2))/2  )        
        }
        
        density.eV <- log(mean(exp(density.eV.holder)))
        if(abs(density.eV) == Inf){
            cat("    Precision reinforced! \n")
            print(density.eV.holder)
            density.eV <- as.numeric(log(mean(exp(mpfr(density.eV.holder, precBits=53)))))
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
            
            ## Step 1: update bhat
            if(constant){
                bhat <- updateb(Z, MU, s2, XtX, b0, B0)
                Zb <- Z - bhat
            }
            
            ## Step 2: update U
            U <- updateU(K, U, V, R, Zb, s2, eU.st, iVU.st)
            if (UL.Normal == "Normal"){
                U <- Unormal(U)
            }else if(UL.Normal == "Orthonormal"){
                U <- GramSchmidt(U)
            }
            
            ## Step 3: update V
            V <- updateV(Zb, U, R, K, s2, eV.st, iVV, UTA)         
            ## update MU
            MU <- M.U(list(U,U,V))
            
            ## Step 4: update eU and iVU
            ## SS <-  diag(psi0, nrow=R) +  t(U)%*%U 
            ## iVU <- rwish(solve(SS), Psi0+K[1])
            ## eU <- c(rMVNorm(1, apply(U,2,sum)/(K[1]+1), solve(iVU)/(K[1]+1)))
            
            ## Step 5: update eV and iVV
            dens.temp <- rep(NA, R)
            SS <-  t(V)%*%V ## (K[3]-1)*cov(V)   +K[3]*msi/(K[3]+1)
            ## invSS <- chol2inv(chol(SS))
            ## iVV <- rwish(invSS, psi.V1)
            ## density.iVV.holder[g] <- dwishartc(U=iVV.st, nu=psi.V1, S=invSS, log=TRUE) ## log(dwish(iVV.st, psi.V1, invSS))                       
            ## eV <- c(rMVNorm(1,apply(V,2,sum)/(K[3]+1), solve(iVV)/(K[3]+1)))
            for(r in 1:R){
                iVV[r,r] <- 1/rgamma(1, (v0 + K[3])/2, (v1+ SS[r,r])/2)
                dens.temp[r] <- log(dinvgamma(iVV.st[r,r], v0/2, v1/2))
            }
            ## invSS <- chol2inv(chol(SS))
            ## iVU[[j]] <- rwish(invSS, psi.U0 + Km[[j]][1])
            density.iVV.holder[g] <- sum(dens.temp)
            
            
            ## Step 6: update s2
            ZE <- Zb - MU 
            s2 <- rinvgamma(1, (c0 + XtX)/2, (d0+ sum(c(ZE)^2))/2  )        
        }
        
        density.iVV <- log(mean(exp(density.iVV.holder)))
        if(abs(density.iVV) == Inf){
            cat("    Precision reinforced! \n")
            print(density.iVV.holder)
            density.iVV <- log(mean(exp(mpfr(density.iVV.holder, precBits=53))))
        }
        
        cat("\n---------------------------------------------- \n ")
        cat("Marignal Likelihood Computation Step 4 \n")
        cat("    density.iVV: ", as.numeric(density.iVV), "\n")
        cat("---------------------------------------------- \n ")

        if(constant){
            ##
            ## Marginal Step 5: p(bhat.st|Z, eU.st, iVU.st, eV.st, iVV.st)
            ##
            density.bhat.holder <- matrix(NA, reduce.mcmc)
            for (g in 1:reduce.mcmc){
                
                ## Step 1: update bhat
                ZU <- Z - MU
                Xtz <- sum(ZU) ## t(apply(X, 4, c))%*%c(ZU)
                cV <- 1/(XtX/s2 +  1/B0) ## solve( XtX/s2 +  diag(1/B0, p))
                cE <- cV*(Xtz/s2 + (1/B0)*b0) ## cV%*%(Xtz/s2 + diag(1/B0, p)*b0)
                bhat <- rnorm(1, cE, sqrt(cV)) ## rMVNorm(1,cE,cV)[1:p]
                Zb <- Z - bhat
                density.bhat.holder[g] <- dnorm(bhat.st, cE, sqrt(cV), log=TRUE)            
                ## bhat <- updateb(Z, MU, s2, XtX, b0, B0)
                ## Zb <- Z - bhat
                
                ## Step 2: update U
                U <- updateU(K, U, V, R, Zb, s2, eU.st, iVU.st)
                if (UL.Normal == "Normal"){
                    U <- Unormal(U)
                }else if(UL.Normal == "Orthonormal"){
                    U <- GramSchmidt(U)
                }
                
                ## Step 3: update V
                V <- updateV(Zb, U, R, K, s2, eV.st, iVV.st, UTA)         
                ## update MU
                MU <- M.U(list(U,U,V))
                
                ## Step 4: update eU and iVU
                ## SS <-  diag(psi0, nrow=R) +  t(U)%*%U 
                ## iVU <- rwish(solve(SS), Psi0+K[1])
                ## eU <- c(rMVNorm(1, apply(U,2,sum)/(K[1]+1), solve(iVU)/(K[1]+1)))
                
                ## Step 5: update eV and iVV
                ## SS <-  diag(psi0, nrow=R) +  t(V)%*%V ## (K[3]-1)*cov(V)   +K[3]*msi/(K[3]+1)
                ## iVV <- rwish(solve(SS), Psi0+K[3])
                ## eV <- c(rMVNorm(1,apply(V,2,sum)/(K[3]+1), solve(iVV)/(K[3]+1)))
                
                ## Step 6: update s2
                ZE <- Zb - MU 
                s2 <- rinvgamma(1, (c0 + XtX)/2, (d0+ sum(c(ZE)^2))/2  )        
            }
            
            density.bhat <- log(mean(exp(density.bhat.holder)))
            if(abs(density.bhat) == Inf){
                cat("    Precision reinforced! \n")
                print(density.bhat.holder)
                density.bhat <- as.numeric(log(mean(exp(mpfr(density.bhat.holder, precBits=53)))))
            }
            cat("\n---------------------------------------------- \n ")
            cat("Marignal Likelihood Computation Step 5 \n")
            cat("    density.bhat: ", as.numeric(density.bhat), "\n")
            cat("---------------------------------------------- \n ")
        }

        ##
        ## Marginal Step 6: p(sigma.st|Z, eU.st, iVU.st, eV.st, iVV.st, bhat.st)
        ##
        density.Sigma.holder <- matrix(NA, reduce.mcmc)
        for (g in 1:reduce.mcmc){
            
            ## Step 1: update bhat
            ## bhat <- updateb(Z, MU, s2, XtX, b0, B0)
            ## Zb <- Z - bhat.st
            
            ## Step 2: update U
            U <- updateU(K, U, V, R, Zb.st, s2, eU.st, iVU.st)
            if (UL.Normal == "Normal"){
                U <- Unormal(U)
            }else if(UL.Normal == "Orthonormal"){
                U <- GramSchmidt(U)
            }
            
            ## Step 3: update V
            V <- updateV(Zb.st, U, R, K, s2, eV.st, iVV.st, UTA)         
            ## update MU
            MU <- M.U(list(U,U,V))
            
            ## Step 4: update eU and iVU
            ## SS <-  diag(psi0, nrow=R) +  t(U)%*%U 
            ## iVU <- rwish(solve(SS), Psi0+K[1])
            ## eU <- c(rMVNorm(1, apply(U,2,sum)/(K[1]+1), solve(iVU)/(K[1]+1)))
            
            ## Step 5: update eV and iVV
            ## SS <-  diag(psi0, nrow=R) +  t(V)%*%V ## (K[3]-1)*cov(V)   +K[3]*msi/(K[3]+1)
            ## iVV <- rwish(solve(SS), Psi0+K[3])
            ## eV <- c(rMVNorm(1,apply(V,2,sum)/(K[3]+1), solve(iVV)/(K[3]+1)))
            
            ## Step 6: update s2
            ZE <- Zb - MU
            nu1 <- (c0 + XtX)/2
            nu2 <- (d0+ sum(c(ZE)^2))/2  
            s2 <- rinvgamma(1, nu1, nu2)
            density.Sigma.holder[g] <- log(dinvgamma(Sigma.st, nu1, nu2))
        }
        
        density.Sigma <- log(mean(exp(density.Sigma.holder)))
        if(abs(density.Sigma) == Inf){
            cat("    Precision reinforced! \n")
            print(density.Sigma.holder)
            density.Sigma <- as.numeric(log(mean(exp(mpfr(density.Sigma.holder, precBits=53)))))
        }
        
        cat("\n---------------------------------------------- \n ")
        cat("Marignal Likelihood Computation Step 6 \n")
        cat("    density.Sigma: ", as.numeric(density.Sigma), "\n")
        cat("---------------------------------------------- \n ")
        
        
        ## Prior ordinate
        iVV <- iVU <- diag(R) ; eV <- eU <- rep(0,R) 
        density.eU.prior <- log(dmvnorm(eU.st, eU, iVU)) ## sum(sapply(1:K[1], function(j){log(dmvnorm(U.st[j,], eU, iVU))}))
        density.eV.prior <- log(dmvnorm(eV.st, eV, iVV)) ## sum(sapply(1:K[3], function(j){log(dmvnorm(L.st[j,], eV, iVV))}))
        ## density.iVU.prior <- dwishartc(U=iVU.st, nu=psi.U0, S=Psi.U0, log=TRUE) ## log(dwish(iVU.st, psi.U0, Psi.U0)) 
        ## density.iVV.prior <- dwishartc(U=iVV.st, nu=psi.V0, S=Psi.V0, log=TRUE) ## log(dwish(iVV.st, psi.V0, Psi.V0))
        dens.temp.iVU <- dens.temp.iVV <- rep(NA, R)
        for(r in 1:R){
            dens.temp.iVU[r] <- log(dinvgamma(iVU.st[r,r], u0/2, u1/2))
            dens.temp.iVV[r] <- log(dinvgamma(iVV.st[r,r], v0/2, v1/2))
        }
        density.iVU.prior <- sum(dens.temp.iVU)
        density.iVV.prior <- sum(dens.temp.iVV)

        if(constant){
            density.beta.prior <- dnorm(bhat.st, b0, sqrt(B0), log=TRUE) ## p = 1
        }
        density.Sigma.prior <- log(dinvgamma(Sigma.st, c0/2, (d0)/2))


        if(constant){
            logprior <- density.eU.prior + density.eV.prior + density.iVU.prior + density.iVV.prior +
                density.beta.prior + density.Sigma.prior
            logdenom <- (density.eU + density.eV + density.iVU + density.iVV +
                         density.bhat + density.Sigma)
        }else{
            logprior <- density.eU.prior + density.eV.prior + density.iVU.prior + density.iVV.prior +
                density.Sigma.prior
            logdenom <- (density.eU + density.eV + density.iVU + density.iVV + density.Sigma)
        }
        ## logmarglike <- (loglike + logprior) - logdenom;
        logmarglike.upper <- (loglike.upper + logprior) - logdenom
        
        cat("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n")
        ## cat("    log marginal likelihood = (loglike + logprior) - (density.parameters) \n")
        ## cat("    log marginal likelihood : ", as.numeric(logmarglike), "\n")
        cat("    log marginal likelihood : ", as.numeric(logmarglike.upper), "\n")
        cat("    logprior: ", as.numeric(logprior), "\n")
        cat("    log posterior density: ", as.numeric(logdenom), "\n")
        cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n")
                
    }
    attr(output, "title") <- "NetworkStatic Posterior Sample"
    attr(output, "Z") <- Z
    attr(output, "m") <- 0 
    attr(output, "eU") <- eUmat
    attr(output, "eV") <- eVmat
    attr(output, "iVU") <- iVUmat
    attr(output, "iVV") <- iVVmat
    attr(output, "U") <- U
    attr(output, "V") <- V
    attr(output, "MU") <- MU.record/nss
    attr(output, "Umat") <- Umat
    attr(output, "Vmat") <- Vmat
    if(constant){
        attr(output, "bmat") <- bmat
    }
    attr(output, "s2mat") <- s2mat
    attr(output, "mcmc") <- nstore
    attr(output, "R") <- R
    attr(output, "DIC") <- Z.DIC
    attr(output, "Waic.out") <- Waic.out
    ## attr(output, "loglike") <- loglike
    attr(output, "loglike") <- loglike.upper
    ## attr(output, "logmarglike") <- logmarglike
    attr(output, "logmarglike") <- logmarglike.upper 
    ## cat("elapsed time is ", proc.time() - ptm, "\n")
    return(output)
}
