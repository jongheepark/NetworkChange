% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MarginalCompare.R
\name{MarginalCompare}
\alias{MarginalCompare}
\title{Compare Log Marginal Likelihood}
\usage{
MarginalCompare(outlist)
}
\arguments{
\item{outlist}{List of NetworkChange objects}
}
\value{
A matrix of log marginal likelihoods.
}
\description{
Compare Log Marginal Likelihood
}
\references{
Siddhartha Chib. 1995. ``Marginal Likelihood from the Gibbs Output.''
\emph{Journal of the American Statistical Association}. 90: 1313-1321.

Jong Hee Park and Yunkyun Sohn. 2020. "Detecting Structural Change
in Longitudinal Network Data." \emph{Bayesian Analysis}. Vol.15, No.1, pp.133-157.

Sumio Watanabe. 2010. "Asymptotic equivalence of Bayes cross validation and widely
applicable information criterion in singular learning theory."
\emph{Journal of Machine Learning Research}. 11: 3571-3594.
}
\seealso{
\code{\link{WaicCompare}}
}
