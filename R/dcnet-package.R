#' The 'dcnet' package.
#'
#' @description Estimate dynamic conditional partial correlation networks for multilevel time-series data. These
#' models are based on multilevel vector autoregressive processes for the means and multilevel dynamic conditional
#' correlations (mlDCC) for the residuals. The mlDCC structure contains the longitudinal partial correlations that
#' are estimated for each time-point. Sampling is done in stan, a C++ package providing HMC methods for full Bayesian
#' inference (cf. [http://mc-stan.org]). either via variational Bayes or Hamiltonian Monte Carlo sampling. 
#'
#' @docType package
#' @name dcnet-package
#' @aliases dcnet
#' @useDynLib dcnet, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#'
NULL
