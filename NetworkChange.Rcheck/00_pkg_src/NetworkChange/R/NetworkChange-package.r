#' @import ggplot2
#' @import graphics
#' @import MASS
#' @import stats
#' @import utils
#' @import grDevices
#' @import methods
#' @import mvtnorm
#' @import RColorBrewer
#' @import MCMCpack
#' @importFrom Rmpfr mpfr
#' @importFrom patchwork wrap_plots
#' @importFrom viridis viridis
NULL
#> NULL

#' Fit Bayesian multilinear tensor regression model with Change-points.
#'
#' This package implements gibbs updates for Bayesian multilinear tensor regression model
#'  with change-points. Version 1.0 includes performance optimizations with vectorized
#'  MCMC operations and modern ggplot2-based visualizations with colorblind-friendly palettes.
#'
