##' Draw samples from a longitudinal partial correlation network with dynamic correlations. 
##'
##' The fitted models are 'cmdstanr' objects and all posterior parameter estimates can be obtained and can be examined with either the 'rstan' toolbox, plotted and printed using generic functions  or passed to 'dcnet' functions to 'forecast' or compute 'model_weights' or compute fit statistics based on leave-future-out cross-validation. 
##' 
##' @title Estimate a Bayesian Dynamic Correlation Network
##' @param data Time-series or matrix object. A time-series or matrix object containing observations at the same interval.
##' @param J ...
##' @param group Vector with group id
##' @param xC Numeric vector or matrix. Covariates(s) for the constant variance terms in C, or c, used in a log-linear model on the constant variance terms \insertCite{Rast2020}{dcnet}. If vector, then it acts as a covariate for all constant variance terms. If matrix, must have columns equal to number of time series, and each column acts as a covariate for the respective time series (e.g., column 1 predicts constant variance for time series 1).
##' @param S_pred Dummy for S, for each time-point (same dimension as data)
##' @param parameterization Character (Default: "DCC"). One of 'CCC', 'DCCr' or 'DCC'. 
##' @param P Integer. Dimension of GARCH component in MGARCH(P,Q).
##' @param Q Integer. Dimension of ARCH component in MGARCH(P,Q).
##' @param iterations Integer (Default: 2000). Number of iterations for each chain (including warmup).
##' @param chains Integer (Default: 4). The number of Markov chains.
##' @param standardize_data Logical (Default: FALSE). Whether data should be standardized. 
##' @param distribution Character (Default: "Student_t"). Distribution of innovation: "Student_t"  or "Gaussian"
##' @param meanstructure Character (Default: "constant"). Defines model for means. Either 'constant', 'VAR', or 'ARMA'. Currently ARMA(1,1) or 'VAR' (VAR1).
##' @param sampling_algorithm Character (Default" "variational"). Define sampling algorithm. Either 'HMC' or variational Bayes 'variational'.
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
                  J,
                  group,
                  xC = NULL,
                  S_pred =  NULL,
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
                  lbound =  FALSE,
                  ubound =  FALSE, ...) {
    if ( tolower(distribution) == "gaussian" ) {
        num_dist <- 0
    } else if ( tolower(distribution) == "student_t" ) {
        num_dist <- 1
    } else {
        stop( "\n\n Specify distribution: Gaussian or Student_t \n\n")
    }

    stan_data <- stan_data(data,
                           J,
                           group,
                           xC,
                           S_pred,
                           P,
                           Q,
                           standardize_data,
                           distribution = num_dist,
                           meanstructure,
                           simplify_ch,
                           simplify_ah,
                           simplify_bh,
                           lbound,
                           ubound)

    ## Select stanmodel
    ## 
    ## stanmodel <- switch(parameterization,
    ##                     CCC = stanmodels$CCCMGARCH,
    ##                     DCC = stanmodels$DCCMGARCH,
    ##                     NULL)

    ## select CmdStanModel created by cmdstan_model at time of compilation
    ## rstan
    stan_path <- .get_target_stan_path()

    ccc_file <- file.path(stan_path, "VAR.stan" )
    dcc_file <- file.path(stan_path, "DCCMGARCHrandQ.stan" )
    dccr_file <-file.path(stan_path, "DCCMGARCHrandS.stan" )
    
    stanmodel <- switch(parameterization,
                        CCC = cmdstan_model(ccc_file, include_paths =  stan_path,
                                            cpp_options = list(stan_threads = TRUE)),
                        DCC = cmdstan_model(dcc_file, include_paths =  stan_path,
                                            cpp_options = list(stan_threads = TRUE)),
                        DCCr = cmdstan_model(dccr_file, include_paths =  stan_path,
                                             cpp_options = list(stan_threads = TRUE)),
                        NULL)
        
    if(is.null(stanmodel)) {
        stop("Not a valid model specification. ",
             parameterization,
             "must be one of: ",
             paste0(supported_models, collapse = ", "),
             ".")
    }

    ## HMC Sampling
    if( tolower( sampling_algorithm ) == 'hmc') {
      if( is.null( iterations ) ) {
        iter_warmup <- 1000
        iter_sampling <- 1000
      } else if ( !is.null( iterations ) ) {
        iter_warmup <- round(iterations/2 )
        iter_sampling <- iterations - iter_warmup
      }

      ## 
      max_cores <- parallel::detectCores()
      model_fit <- stanmodel$sample(data = stan_data,
                                    iter_warmup = iter_warmup,
                                    iter_sampling =  iter_sampling,
                                    adapt_delta = .95,
                                    chains = chains,
                                    parallel_chains = min( max_cores,  chains ), 
                                    ...)
    } else if ( tolower(sampling_algorithm) == 'variational' ) {
      ## Sampling via Variational Bayes
      if( is.null( iterations ) ) iterations <- 30000
      model_fit <- stanmodel$variational(data = stan_data,
                                         iter = iterations, ...)
    } else {
      stop( "\n\n Provide sampling algorithm: 'HMC' or 'variational'\n\n" )
    }
    
    ## Model fit is based on standardized values.
    if(standardize_data ) {
      mns <- stan_data$grand_mean
      sds <- stan_data$grand_sd
    } else {
      mns <- stan_data$grand_mean
      sds <- 1
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
                       sampling_algorithm = sampling_algorithm)

    class(return_fit) <- "dcnet"
    return(return_fit)
}


#' Models supported by dcnet
#'
#' To be used when checking whether a parameterization or object type is a supported type.
#' May facilitate more parameterizations, as we only have to update these, and the switch statements.
#' @keywords internal
#' @author Philippe Rast
supported_models <- c("CCC", "DCC", "DCCr")
