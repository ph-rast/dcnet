% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dcnet.R
\name{dcnet}
\alias{dcnet}
\title{Estimate a Bayesian Dynamic Correlation Network}
\usage{
dcnet(
  data,
  J = NULL,
  group = NULL,
  xC = NULL,
  S_pred = NULL,
  parameterization = "DCC",
  P = 1,
  Q = 1,
  iterations = NULL,
  chains = 4,
  standardize_data = FALSE,
  distribution = "Gaussian",
  meanstructure = "VAR",
  sampling_algorithm = "variational",
  simplify_ch = 1,
  simplify_ah = 1,
  simplify_bh = 1,
  lbound = FALSE,
  ubound = FALSE,
  threads_per_chain = 4,
  grainsize = 1,
  ...
)
}
\arguments{
\item{data}{Time-series or matrix object. A time-series or matrix object containing observations at the same interval.}

\item{J}{Number of clusters. If \code{NULL}, derived from data, otherwise user provided.}

\item{group}{Vector with group id. If \code{NULL}, derived from data, otherwise user provided.}

\item{xC}{Numeric vector or matrix. Covariates(s) for the constant variance terms in C, or c, used in a log-linear model on the constant variance terms \insertCite{Rast2020}{dcnet}. If vector, then it acts as a covariate for all constant variance terms. If matrix, must have columns equal to number of time series, and each column acts as a covariate for the respective time series (e.g., column 1 predicts constant variance for time series 1).}

\item{S_pred}{Dummy for S, for each time-point (same dimension as data)}

\item{parameterization}{Character (Default: "DCC"). One of 'CCC', 'DCCr' or 'DCC', or 'DCCrs' for multithread.}

\item{P}{Integer. Dimension of GARCH component in MGARCH(P,Q).}

\item{Q}{Integer. Dimension of ARCH component in MGARCH(P,Q).}

\item{iterations}{Integer (Default: 2000). Number of iterations for each chain (including warmup).}

\item{chains}{Integer (Default: 4). The number of Markov chains.}

\item{standardize_data}{Logical (Default: FALSE). Whether data should be standardized.}

\item{distribution}{Character (Default: "Student_t"). Distribution of innovation: "Student_t"  or "Gaussian"}

\item{meanstructure}{Character (Default: "constant"). Defines model for means. Either 'constant', 'VAR', or 'ARMA'. Currently ARMA(1,1) or 'VAR' (VAR1).}

\item{sampling_algorithm}{Character (Default" "variational"). Define sampling algorithm: One of Hamilton Monte-Carlo 'HMC', variational Bayes 'variational' or 'pathfinder'.}

\item{simplify_ch}{Random efx on ch}

\item{simplify_ah}{Randon efx on ah}

\item{simplify_bh}{Random efx on bh}

\item{lbound}{Lower bound on observed values (across data matrix) - defaults to min(data).}

\item{ubound}{Upper bound on observed values (across data matrix) - defaults to max(data).}

\item{...}{Additional arguments can be ‘chain_id’, ‘init_r’, ‘test_grad’, ‘append_samples’, ‘refresh’, ‘enable_random_init’ etc.}
}
\value{
\code{dcnet} object.
}
\description{
Draw samples from a longitudinal partial correlation network with dynamic correlations.
}
\details{
The fitted models are 'cmdstanr' objects and all posterior parameter estimates can be obtained and can be examined with either the 'rstan' toolbox, plotted and printed using generic functions  or passed to 'dcnet' functions to 'forecast' or compute 'model_weights' or compute fit statistics based on leave-future-out cross-validation.
}
\examples{
\dontrun{
##
}
}
\references{
\insertAllCited()
}
\author{
Philippe Rast
}
