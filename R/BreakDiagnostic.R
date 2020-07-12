#' Detect a break number using different metrics 
#'
#'
#' @param Y Reponse tensor 
#' @param R Dimension of latent space. The default is 2. 
#' @param break.upper Upper threshold for break number detection.
#'  The default is \code{break.upper = 3}.
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
#'
#' @param UL.Normal Transformation of sampled U. Users can choose "NULL", "Normal" or "Orthonormal."
#' "NULL" is no normalization. "Normal" is the standard normalization.
#' "Orthonormal" is the Gram-Schmidt orthgonalization. Default is "NULL."
#'
#' 
#' @param v0 \eqn{v_0/2} is the shape parameter for the inverse
#' Gamma prior on variance parameters for V.
#' If \code{v0 = NULL}, a value is computed from a test run of \code{NetworkStatic}.
#' 
#' @param v1 \eqn{v_1/2} is the scale parameter for the
#' inverse Gamma prior on variance parameters for V.
#' If \code{v1 = NULL}, a value is computed from a test run of \code{NetworkStatic}.
#'
#' @param a \eqn{a} is the shape1 beta prior for transition probabilities. By default,
#' the expected duration is computed and corresponding a and b values are assigned. The expected
#' duration is the sample period divided by the number of states.
    
#' @param b \eqn{b} is the shape2 beta prior for transition probabilities. By default,
#' the expected duration is computed and corresponding a and b values are assigned. The expected
#' duration is the sample period divided by the number of states.
#'
#' @importFrom tidyr gather
#'
#' @references   Jong Hee Park and Yunkyun Sohn. 2020. "Detecting Structural Change
#' in Longitudinal Network Data." \emph{Bayesian Analysis}. Vol.15, No.1, pp.133-157.
#' @export
#'
#'
#' @examples
#'    \dontrun{
#'    set.seed(19333)
#'    ## Generate an array (15 by 15 by 20) with a block merging transition
#'    Y <- MakeBlockNetworkChange(n=5, T=20, type ="merge")
#' 
#'    ## Fit 3 models (no break, one break, and two break) for break number detection 
#'    detect <- BreakDiagnostic(Y, R=2, break.upper = 2)
#'    
#'    ## Look at the graph
#'    detect[[1]]; print(detect[[2]])
#'
#' }
#'
#'
#'

BreakDiagnostic <- function(Y, R=2, mcmc=100, burnin=100, verbose=100, thin=1, UL.Normal = "Orthonormal",
                            v0=NULL, v1=NULL, break.upper = 3, a=1, b=1){
    ## set.seed(11173)
    K <- dim(Y)
    ## prior estimate
    if(is.null(v0) || is.null(v1)){
        test.run <- NetworkStatic(Y=Y, R=R, mcmc=mcmc, burnin=burnin, thin=thin, verbose=verbose, UL.Normal = UL.Normal,
                                  v0=10, v1=K[[3]]*2)
        V <- attr(test.run, "V")
        sigma.mu = mean(apply(V, 2, mean))
        sigma.var = 10*mean(apply(V, 2, var))
        if(is.null(v0)) v0 <- 4 + 2 * (sigma.mu^2/sigma.var)
        if(is.null(v1)) v1 <- 2 * sigma.mu * (v0/2 - 1)
    }
    
    ## model fit
    out <- as.list(rep(NA, break.upper))
    out[[1]] <- NetworkStatic(Y, R=2, mcmc=mcmc, reduce.mcmc = mcmc/2, 
                              burnin=burnin, verbose=verbose, v0=v0, v1=v1,
                              Waic=TRUE, marginal=TRUE)
    for(m in 1:break.upper){
        ## to save time and to be more conservative, use randomly generated initial states
        initial.s <- sort(rep(1:(m+1), length=K[[3]]))
        out[[m+1]] <- NetworkChange(Y, R=2, m=m, mcmc=mcmc, initial.s = initial.s, reduce.mcmc = mcmc/2, 
                                    burnin=burnin, verbose=verbose, thin=thin, v0=v0, v1=v1, a=a, b=a,
                                    Waic=TRUE, marginal=TRUE)
    }
    
    ## diagnostic info
    Waic.holder <- marginal.holder <- loglike.holder <- rep(NA, break.upper+1)
    for(i in 1:(break.upper+1)){
        loglike.holder[i] <- attr(out[[i]], "loglike")
        marginal.holder[i] <- attr(out[[i]], "logmarglike")
        Waic.holder[i] <- attr(out[[i]], "Waic.out")[1]
    }
    ## loss
    loss.input <- out[-1]
    loss.out <- BreakPointLoss(loss.input, display=FALSE)[[1]]

    ## save model diagnostics
    result <- list("LogMarginal" = marginal.holder,
                   "Loglike" = loglike.holder,
                   "WAIC" = Waic.holder,
                   "Average Loss" = loss.out)

    test.curve1 <- -2*matrix(result[[1]], 1, break.upper +1, byrow=TRUE)
    test.curve2 <- -2*matrix(result[[2]], 1, break.upper +1, byrow=TRUE)
    test.curve3 <- matrix(result[[3]], 1, break.upper +1, byrow=TRUE)
    test.curve4 <- matrix(c(NA, result[[4]]), 1, break.upper +1, byrow=TRUE)

    test.curve <- rbind(test.curve1, test.curve2, test.curve3, test.curve4)
    test.curve <- data.frame(test.curve)
    colnames(test.curve) <- paste0("break", 0:break.upper)
    Metric <- c("-2*LogMarginal", "-2*Loglike", "WAIC", "Average Loss")
    test.curve$Metric <- Metric
    data_long <- tidyr::gather(test.curve, model, value, 
                               paste0("break", 0:break.upper), factor_key=TRUE)
    model <- value <- NULL
    g1 <- ggplot(data= transform(data_long,
                                 Metric=factor(Metric,levels=c("-2*LogMarginal", "-2*Loglike", "WAIC", "Average Loss"))),
                 mapping = aes(x = model, y = value, group = Metric, color = Metric)) +
        geom_line(size=0.2) + geom_point(cex=3, alpha=1/2) + facet_wrap(~Metric, nrow=1, ncol=4, scales = "free_y") + 
        labs(x = "Model", y = "Value") + theme_bw() +
        theme(legend.position="none",
              legend.key = element_blank(),
              plot.title = element_text(hjust = 0.5))  
    return(list(graph=g1, result=result))
}
