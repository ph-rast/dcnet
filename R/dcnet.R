##' Draw samples from a longitudinal partial correlation network with dynamic correlations. 
##'
##' The fitted models are 'cmdstanr' objects and all posterior parameter estimates can be obtained and can be examined with either the 'rstan' toolbox, plotted and printed using generic functions  or passed to 'dcnet' functions to 'forecast' or compute 'model_weights' or compute fit statistics based on leave-future-out cross-validation.
##'
##' The two stage model is the default, where the uncondidional coorelation S and its random effects are estimated from the data and passed to the mlVAR-DCC model.
##' 
##' @title Estimate a Bayesian Dynamic Correlation Network
##' @param data Time-series or matrix object. A time-series or matrix object containing observations at the same interval.
##' @param J Number of clusters. If `NULL`, derived from data, otherwise user provided. 
##' @param group Vector with group id. If `NULL`, derived from data, otherwise user provided. 
##' @param xC Numeric vector or matrix. Covariates(s) for the constant variance terms in C, or c, used in a log-linear model on the constant variance terms \insertCite{Rast2020}{dcnet}. If vector, then it acts as a covariate for all constant variance terms. If matrix, must have columns equal to number of time series, and each column acts as a covariate for the respective time series (e.g., column 1 predicts constant variance for time series 1).
##' @param S_pred Dummy for S, for each time-point (same dimension as data)
##' @param parameterization Character (Default: "DCC"). One of 'CCC', 'DCCr' or 'DCC', or 'DCCrs' for multithread.
##' @param twostage Logical (Defaults to TRUE). Should a two stage approach be used. 
##' @param P Integer. Dimension of GARCH component in MGARCH(P,Q).
##' @param Q Integer. Dimension of ARCH component in MGARCH(P,Q).
##' @param iterations Integer (Default: 2000). Number of iterations for each chain (including warmup).
##' @param chains Integer (Default: 4). The number of Markov chains.
##' @param standardize_data Logical (Default: FALSE). Whether data should be standardized. 
##' @param distribution Character (Default: "Student_t"). Distribution of innovation: "Student_t"  or "Gaussian"
##' @param meanstructure Character (Default: "constant"). Defines model for means. Either 'constant', 'VAR', or 'ARMA'. Currently ARMA(1,1) or 'VAR' (VAR1).
##' @param sampling_algorithm Character (Default" "variational"). Define sampling algorithm: One of Hamilton Monte-Carlo 'HMC', variational Bayes 'variational' or 'pathfinder'.
##' @param simplify_ch Random efx on ch 
##' @param simplify_ah Randon efx on ah
##' @param simplify_bh Random efx on bh
##' @param lbound Lower bound on observed values (across data matrix) - defaults to min(data).
##' @param ubound Upper bound on observed values (across data matrix) - defaults to max(data).
##' @param ... Additional arguments can be ‘chain_id’, ‘init_r’, ‘test_grad’, ‘append_samples’, ‘refresh’, ‘enable_random_init’ etc. 
##' @return \code{dcnet} object.
##' @author Philippe Rast
##' @references
##'    \insertAllCited()
##' @export
##' @examples
##' \dontrun{
##' ##
##' }
dcnet <- function(data,
                  J = NULL,
                  group = NULL,
                  xC = NULL,
                  S_pred =  NULL,
                  parameterization = "DCC",
                  twostage = TRUE,
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
                  lbound =  FALSE,
                  ubound =  FALSE,
                  threads_per_chain = 4,
                  grainsize = 1, ...) {
  
  ## Identify distribution type
  num_dist <- switch(tolower(distribution),
                     "gaussian" = 0,
                     "student_t" = 1,
                     stop("Specify distribution: Gaussian or Student_t"))

  ## Set number of participants if not provided
  if(is.null(J)) J <- N <- length(data)

  ## Create group if not provided
  if(is.null(group)) group <- rep(seq_len(N),  each = nrow(data[[1]]))

  ## Prepare Stan data
  stan_data <- stan_data(data, J, group, xC, S_pred, P, Q,
                         standardize_data, distribution = num_dist,
                         meanstructure, simplify_ch, simplify_ah,
                         simplify_bh, lbound, ubound, grainsize)

    ## Select the correct pre-compiled model from the global environment
    stanmodel <- switch(parameterization,
                        CCC = ccc_model,
                        DCC = dcc_model,
                        DCCr = dccr_model,
                        DCCrs = dccrs_model,
                        stop("Invalid parameterization"))
    
    ## Precompute S if twostage == TRUE
    if(twostage == TRUE) {
        ## Compute individual sample correlations
        ## Extract fisher z transform off diagional elements
        zhat <- lapply(seq_len(J), function(x) {
            Rj <- cor(data[[x]])
            atanh(Rj[lower.tri(Rj)])
        })
        ## S data
        nt <- stan_data$nt
        Sdim <-  nt*(nt-1)/2
        s_data <- list(J = J,
                       nt = stan_data$nt,
                       Sdim = Sdim,
                       zhat = zhat) 
        if (is.null(iterations)) iterations <- 30000
        cat("Stage 1: Random effect of S is estimated \n")
        precomp_fit <- stage2_model$variational(data = s_data,
                                                iter = iterations,
                                                threads = threads_per_chain,
                                                ...)
        cat("\n Stage 2: mlVAR-DCC is estimated \n")
        ## Extract the random effects SD vector
        post2draws <- precomp_fit$draws(format = "df")
        sigma_z <- post2draws$sigma_z
        ## Simplest approach: Average across draws and use E(sigma_z) as SD for all S ranefs
        sigma_z_hat <- rep(mean(sigma_z), Sdim)
        ## TODO expand and include full posterior of sigma_z instead of point estimate
        sigma_z_hat_full <- sigma_z
        ## Attach to stan_data
        stan_data$Sdim <- Sdim
        stan_data$S_vec_tau_fixed <- sigma_z_hat
        stan_data$S_vec_tau_full <- sigma_z_hat_full ## save draw to report posterior Cri etc. 
    }
    
  ## if(is.null(stanmodel)) {
  ##     stop("Not a valid model specification. ",
  ##          parameterization,
  ##          "must be one of: ",
  ##          paste0(supported_models, collapse = ", "),
  ##          ".")
  ## }

    ## HMC Sampling
    if (tolower(sampling_algorithm) == "hmc") {
        if (is.null(iterations)) {
            iter_warmup <- 1000
            iter_sampling <- 1000
        } else if (!is.null(iterations)) {
            iter_warmup <- round(iterations / 2)
            iter_sampling <- iterations - iter_warmup
        }

        ##
        max_cores <- parallel::detectCores()
        Sys.setenv(STAN_NUM_THREADS = threads_per_chain)
        model_fit <- stanmodel$sample(data = stan_data,
                                      iter_warmup = iter_warmup,
                                      iter_sampling = iter_sampling,
                                      adapt_delta = .95,
                                      chains = chains,
                                      parallel_chains = min(max_cores, chains),
                                      threads_per_chain = threads_per_chain,
                                      ...)
    } else if (tolower(sampling_algorithm) == "variational") {
        ## Sampling via Variational Bayes
        if (is.null(iterations)) iterations <- 30000
        model_fit <- stanmodel$variational(data = stan_data,
                                           iter = iterations,
                                           threads = threads_per_chain,
                                          ...)
    } else if (tolower(sampling_algorithm) == "pathfinder") {
        max_cores <- parallel::detectCores()
        Sys.setenv(STAN_NUM_THREADS = threads_per_chain)
        ## Sampling via Pathfinder Method
        if (is.null(iterations)) iterations <- 30000
        model_fit <- stanmodel$pathfinder(data = stan_data,
                                          num_threads = threads_per_chain,
                                        #iter = iterations,
                                        #jacobian =  TRUE,
                                          ...)
    } else if (tolower(sampling_algorithm) == "laplace") {
        max_cores <- parallel::detectCores()
        Sys.setenv(STAN_NUM_THREADS = threads_per_chain)
        fit_mode <- stanmodel$optimize(data = stan_data, jacobian = TRUE, threads = threads_per_chain, ...)
        model_fit <- stanmodel$laplace(data = stan_data, mode = fit_mode, threads = threads_per_chain, ...)
    } else {
        stop("\n\n Provide sampling algorithm: 'HMC', 'variational' or 'pathfinder' \n\n")
    }

  ## Model fit is based on standardized values.
  if (standardize_data) {
    mns <- stan_data$grand_mean
    sds <- stan_data$grand_sd
  } else {
    mns <- stan_data$grand_mean
    sds <- 1
  }

  ## Pass out information on whether S_pred is not null. This is for print.R
  if (!is.null(S_pred)) {
    S_pred <- "present"
  }

  ## Values could be converted to original scale using something like this on the estimates
  ## orig_sd = stan_data$rts %*% diag(sds)
  ## orig_scale = orig_sd + array(rep(mns, each = aussi[[1]]$T), dim = c(aussi[[1]]$T, aussi[[1]]$nt) )
  return_fit <- list(model_fit = model_fit,
                     param = parameterization,
                     distribution = distribution,
                     num_dist = num_dist,
                     iter = iterations,
                     chains = chains,
                     elapsed_time = model_fit$time()$total,
                     date = date(),
                     nt = stan_data$nt,
                     TS_length = stan_data$T,
                     TS_names = colnames(stan_data$rts[[1]]),
                     ##RTS_last = stan_data$rts[stan_data$T,],
                     grand_mean = mns,
                     grand_sd = sds,
                     RTS_full = stan_data$rts,
                     mgarchQ = stan_data$Q,
                     mgarchP = stan_data$P,
                     xC = stan_data$xC,
                     meanstructure = stan_data$meanstructure,
                     std_data = standardize_data,
                     sampling_algorithm = sampling_algorithm,
                     S_pred = S_pred)

  class(return_fit) <- "dcnet"
  return(return_fit)
}


#' Models supported by dcnet
#'
#' To be used when checking whether a parameterization or object type is a supported type.
#' May facilitate more parameterizations, as we only have to update these, and the switch statements.
#' @keywords internal
#' @author Philippe Rast
supported_models <- c("CCC", "DCC", "DCCr", "DCCrs")
