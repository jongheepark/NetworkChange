% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/updateV.R
\name{updateV}
\alias{updateV}
\title{Update layer specific network generation rules}
\usage{
updateV(Zb, U, R, K, s2, eV, iVV, UTA)
}
\arguments{
\item{Zb}{Z - beta.}

\item{U}{The latent node positions.}

\item{R}{The dimension of latent space.}

\item{K}{The dimension of Z.}

\item{s2}{The variance of error.}

\item{eV}{The mean of V.}

\item{iVV}{The variance of V.}

\item{UTA}{Indicator of upper triangular array}
}
\value{
A matrix of layer specific network generation rules
}
\description{
Update layer specific network generation rules
}
