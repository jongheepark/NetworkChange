#' k-fold Cross-Validation for Network Change point models
#'
#' NetworkChangeCVtest implements a suite of k-fold cross-validation test for network change point models

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
#' @param k.fold The number of folds that should be used for cross-validation.
#' @param LOO If \code{TRUE}, do the leave-one-out cross-validation.
#' @param LOO.bound The boundary of LOO test. By default, \code{LOO.bound=2}, meaning that first and last 2 observations are dropped
#' in the leave-one-out cross-validation.
#' @param break.upper The upper limit of break number for cross-validation test.
#' 
#' @return Mean squared error of the k-fold cross validation test 
#'
#' @seealso \code{\link{NetworkChange}}
#'
#' @references  Jong Hee Park and Yunkyun Sohn. 2017. "Detecting Structural Change
#' in Network Time Series Data using Bayesian Inference." Working Paper.
#'
#' @export
#' @examples
#'
#'    \dontrun{
#'    set.seed(1973)
#'
#'    ## generate an array with three constant blocks
#'    set.seed(11173)
#'    N <- 10
#'    Y <- MakeBlockNetworkChange(n=N, break.point = .5, base.prob=.2, block.prob=.5,
#'         shape=1, T=20, type ="split")

#'    G <- 10 ## Small mcmc scans to save time
#'    MSE <- NetworkChangeCVtest(Y, k.fold = 2)
#'
#'    ## find mean of MSE
#'    apply(MSE, 2, mean)
#' 
#'    }
#'

NetworkChangeCVtest <- function(Y, R=2, k.fold = 2, mcmc = 100, burnin=100, thin=1, verbose=0,
                                LOO=FALSE, LOO.bound=2, break.upper = 3, Waic=FALSE){

    K <- dim(Y)
    ntime <- K[3]
    bound <- LOO.bound
    
    ## leave one out test
    if(LOO){
        progress.bar <- create_progress_bar("text")
        test.length <- length(bound:(ntime-bound))
        progress.bar$init(test.length)
        MSE <- matrix(NA, test.length, break.upper+1)
        count <- 1
        for(t in bound:(ntime-bound)){
            y.train <- Y[,,-t]
            y.test <- Y[,,t]
            Time <- dim(y.train)[3]

            ## list object storage
            m.out <- mse <- as.list(rep(NA, break.upper))

            ## fit a break specific model
            m0 <- NetworkStatic(y.train, R=R,  mcmc=mcmc, burnin=burnin, verbose=verbose)
            for(m in 1:break.upper){
                m.out[[m]] <- NetworkChange(y.train, m = m, R=R,  mcmc=mcmc, burnin=burnin, verbose=verbose,
                                            initial.s = sort(rep(1:(m+1), length=Time)))
            }

            ## compute z.test
            Z.test <- y.test
            ee <- eigen(y.test)
            Z.test <- y.test - ee$values[1] * outer(ee$vectors[,1], ee$vectors[,1])

            ## compute mse
            mse0 <- sum((Z.test - apply(m0, 1:2, mean))^2)
            for(m in 1:break.upper){
                mse[[m]] <- sum((Z.test - apply(m.out[[m]], 1:2, mean))^2)
            }

            ## save
            MSE[count, ] <- c(mse0, unlist(mse))
            count <- count + 1
            progress.bar$step()
        }
        
    }else{ ## k-fold CV
        progress.bar <- create_progress_bar("text")
        progress.bar$init(k.fold)
        train.ratio <- 1/k.fold
        MSE <- matrix(NA, k.fold, break.upper+1)
        count <- 1
        for(i in 1:k.fold){
            ## train set
            train_rows <- sort(sample(1:ntime, train.ratio*ntime, prob=rep(1/ntime, ntime)))
            y.train <- Y[,,train_rows]
            Time <- dim(y.train)[3]
            
            ## test set
            y.test <- Y[,,-train_rows]
            Z.test <- y.test
            K.test <- dim(y.test)
            for(k in 1:K.test[3]){
                ee <- eigen(y.test[,,k])
                Z.test[,,k] <- y.test[,,k] - ee$values[1] * outer(ee$vectors[,1], ee$vectors[,1])
            }
            
            ## list object storage
            m.out <- mse <- as.list(rep(NA, break.upper))

            ## fit a break specific model
            m0 <- NetworkStatic(y.train, R=R, mcmc=mcmc, burnin=burnin, verbose=verbose)
            for(m in 1:break.upper){
                m.out[[m]] <- NetworkChange(y.train, m = m, R=R,  mcmc=mcmc, burnin=burnin, verbose=verbose,
                                            initial.s = sort(rep(1:(m+1), length=Time)))
            }
            ## compute mse
            mse0 <- sum((apply(Z.test, 1:2, mean) - apply(m0, 1:2, mean))^2)
            for(m in 1:break.upper){
                mse[[m]] <- sum((apply(Z.test, 1:2, mean) - apply(m.out[[m]], 1:2, mean))^2)
            }

            ## save
            MSE[count, ] <- c(mse0, unlist(mse))
            count <- count + 1
            progress.bar$step()
        }
    }
    if(Waic){
        WAIC.out <- WAIC.list <- as.list(rep(NA, break.upper))
        ## fit a break specific model
        WAIC0 <- NetworkStatic(Y, R=R,  mcmc=mcmc, burnin=burnin, verbose=verbose, Waic=TRUE)
        for(m in 1:break.upper){
            WAIC.out[[m]] <- NetworkChange(Y, m = m, R=R,  mcmc=mcmc, burnin=burnin, verbose=verbose,
                                           initial.s = sort(rep(1:(m+1), length=Time)), Waic=TRUE)
        }
        
        waic0 <- attr(WAIC0, "Waic.out")[1]
        for(m in 1:break.upper){
            WAIC.list[[m]] <- attr(WAIC.out[[m]], "Waic.out")[1]
        }
        
        attr(MSE, "Waic") <- c(waic0, unlist(WAIC.list))
    }
    return(MSE)
}
