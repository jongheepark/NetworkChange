## ---------------------------------------------------- ##
## ---------------------------------------------------- ##
## The idea is to compute the average loss of sampled states from the expected break points
## by JHP "Fri Jul 20 09:57:45 2018"
## ---------------------------------------------------- ##
## ---------------------------------------------------- ##

mse <- function(x, mu){mean((x - mu)^2)}

trace.break <- function(tau.samp, tau){
  m <- length(tau)
  ## par(mfrow=c(1, m))
  if(m>1){
    ## for(i in 1:m){
    ##   plot(tau.samp[i,], type="l"); abline(h=tau[i], col="red")
    ## }
    out <- sapply(1:m, function(i){mse(tau.samp[i,], tau[i])})
  }else{
    ## plot(tau.samp, type="l"); abline(h=tau, col="red")
    out <- mse(tau.samp, tau)
  }
  return(out)
}


findBreakPoint <- function (mcmcout, start = 1) 
{
    out <- attr(mcmcout, "prob.state")
    ## y <- attr(mcmcout, "y")
    m <- attr(mcmcout, "m")
    out[1,] <- c(1, rep(0, m))
   ## if (!is.ts(y)) 
    ##     y <- ts(y, start)
    ## time.frame <- as.vector(time(y))
    if (m == 1) {
        pr.st <- c(0, diff(out[, (m + 1)]))
        pr.st[pr.st < 0] <- 0
        cp <- which(cumsum(pr.st) > 0.5)[1] - 1
    }else {
        cp <- rep(NA, m)
        for (i in 2:m) {
            pr.st <- c(0, diff(out[, i]))
            pr.st <- ifelse(pr.st < 0, 0, pr.st)
            cp[i - 1] <- which(cumsum(pr.st) > 0.5)[1] - 1
        }
        pr.st <- c(0, diff(out[, (m + 1)]))
        pr.st[pr.st < 0] <- 0
        cp[m] <- which(cumsum(pr.st) > 0.5)[1] - 1
    }
    if(sum(is.na(cp))>0){
        cat("\n At break = ", m, " one state is dominated by other states and a break point is not defined for this state. \n")
    }
 
    ## cp.means <- rep(NA, m + 1)
    ## cp.start <- c(1, cp + 1)
    ## cp.end <- c(cp, length(y))
    return(cp + 1)
}


#' Compute the Average Loss of Hidden State Changes from Expected Break Points
#'
#'
#' @param ... MCMC output objects. These have to be of class
#'   \code{mcmc} and have a \code{logmarglike} attribute. In what
#'   follows, we let \code{M} denote the total number of models to be
#'   compared.
#'
#' @param marginal If \code{marginal} is TRUE, \code{logmarglike} will be reported.
#'
#' @param display If \code{display} is TRUE, a plot of \code{ave.loss} will be produced. 
#' 
#' \code{BreakPointLoss}. ave.loss, logmarglike, State, Tau, Tau.samp
#' @return \code{BreakPointLoss} returns five objects. They are: \code{ave.loss} the expected loss for each model
#'   computed by the mean sqaured distance of hidden state changes from the expected break points
#' \textrm{Average Loss} = \frac{1}{M}\sum_{m=1}^{M}\left(\frac{1}{G}\sum_{g=1}^{G} (\bar{\tau}_m - \tau_{m}^{(g)})^2 \right);
#'   \code{logmarglike} the natural log of the marginal likelihood for each model; \code{State} sampled state vectors;
#'   \code{Tau} expected break points for each model; and \code{Tau.samp} sampled break points from hidden state draws.
#'
#' @export
#'
#'
#' @examples
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
#'
#'    ## The most probable model given break number 0 to 3 and data is out1 according to WAIC 
#'    out <- BreakPointLoss(out0, out1, out2, out3, waic=TRUE)
#'
#'    print(out[["ave.loss"]])
#' }
#'
#' 
BreakPointLoss <- function(model.list, waic=FALSE, display=TRUE){

    ## model.list <- list(...)
    M <- length(model.list)
    ## this.call <- match.call()
    ## this.call.string <- deparse(this.call)
    ## this.call.string <- strsplit(this.call.string, "BreakPointLoss\\(")
    ## this.call.string <- this.call.string[[1]][length(this.call.string[[1]])]
    ## this.call.string <- strsplit(this.call.string, ",")

    break.number <- rep(NA, M)
    for (i in 1:M) {
        break.number[i] <- attr(model.list[[i]], "m")
        ## print(break.number[i])
        if(break.number[i] < 1){
            stop("no break model must be dropped\n")
        }
    }

    model.names <- paste0("break ", break.number)## c(model.names, this.call.string[[1]][i])

    ## If marginal, report marginal likelihood
    WAIC <-  NULL
    if(waic){
        WAIC <- WaicCompare(model.list)
    }
    State <- Tau <- Tau.samp <- as.list(rep(NA, M))

    for (i in 1:M) {
        State[[i]] <- attr(model.list[[i]], "Smat")
        Tau[[i]] <- findBreakPoint(model.list[[i]])
        Tau.samp[[i]] <- sapply(1:nrow(State[[i]]), function(j){which(diff(State[[i]][j,])==1)+1})
    }
    
    ## Report Average Loss
    ave.loss <- rep(NA, M)
    for (i in 1:M) {
        ave.loss[i] <- mean(trace.break(Tau.samp[[i]], Tau[[i]]))
    }

    if(display){
        plot(ave.loss, type="o", xlab="Model", ylab="Loss", main="Average Loss", 
             axes=FALSE)
        axis(1, at=break.number, labels=model.names); axis(2)
        abline(v = which.min(ave.loss), lty=3, col="red")
    }
    return(list(ave.loss = ave.loss, WAIC=WAIC, State=State, Tau=Tau, Tau.samp=Tau.samp))
}
