% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BreakPointLoss.R
\name{BreakPointLoss}
\alias{BreakPointLoss}
\title{Compute the Average Loss of Hidden State Changes from Expected Break Points}
\usage{
BreakPointLoss(model.list, waic = FALSE, display = TRUE)
}
\arguments{
\item{model.list}{MCMC output objects. These have to be of class
\code{mcmc} and have a \code{logmarglike} attribute. In what
follows, we let \code{M} denote the total number of models to be
compared.}

\item{waic}{If \code{waic} is TRUE, \code{waic}(Watanabe information criterion) will be reported.}

\item{display}{If \code{display} is TRUE, a plot of \code{ave.loss} will be produced.

\code{BreakPointLoss}. ave.loss, logmarglike, State, Tau, Tau.samp}
}
\value{
\code{BreakPointLoss} returns five objects. They are: \code{ave.loss} the expected loss for each model
computed by the mean sqaured distance of hidden state changes from the expected break points.
\code{logmarglike} the natural log of the marginal likelihood for each model; \code{State} sampled state vectors;
\code{Tau} expected break points for each model; and \code{Tau.samp} sampled break points from hidden state draws.
}
\description{
Compute the Average Loss of Hidden State Changes from Expected Break Points
}
\examples{
   \dontrun{
   set.seed(1973)
   ## Generate an array (30 by 30 by 40) with block transitions
   from 2 blocks to 3 blocks
   Y <- MakeBlockNetworkChange(n=10, T=40, type ="split")
   G <- 100 ## Small mcmc scans to save time

   ## Fit multiple models for break number detection using Bayesian model comparison
   out0 <- NetworkStatic(Y, R=2, mcmc=G, burnin=G, verbose=G, Waic=TRUE)
   out1 <- NetworkChange(Y, R=2, m=1, mcmc=G, burnin=G, verbose=G, Waic=TRUE)
   out2 <- NetworkChange(Y, R=2, m=2, mcmc=G, burnin=G, verbose=G, Waic=TRUE)
   out3 <- NetworkChange(Y, R=2, m=3, mcmc=G, burnin=G, verbose=G, Waic=TRUE)

   ## The most probable model given break number 0 to 3 and data is out1 according to WAIC 
   out <- BreakPointLoss(out0, out1, out2, out3, waic=TRUE)

   print(out[["ave.loss"]])
}


}
\references{
Jong Hee Park and Yunkyun Sohn. 2020. "Detecting Structural Change
in Longitudinal Network Data." \emph{Bayesian Analysis}. Vol.15, No.1, pp.133-157.
}
