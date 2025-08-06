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
##' @param parameterization Character (Default: "DCC"). One of 'CCC', 'DCCj' or 'DCC', or 'DCCms' for multithread.
##' @param multistage Logical (Defaults to TRUE). Should a multi stage approach be used.
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
                  multistage = FALSE,
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

    ## Force into multistage, when DCCms is selected
    if(parameterization == 'DCCms') multistage <- TRUE
    
    ## Identify distribution type
    num_dist <- switch(tolower(distribution),
        "gaussian" = 0,
        "student_t" = 1,
        stop("Specify distribution: Gaussian or Student_t")
    )

    ## Set number of participants if not provided
    if (is.null(J)) J <- N <- length(data)

    ## Create group if not provided
    if (is.null(group)) group <- rep(seq_len(N), each = nrow(data[[1]]))

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
                        DCCms = dccms_model, ## Multistage TODO: Link to multistage argument
                        DCCj = dccj_model,   ## joint TODO: Link to multistage argument
                        stop("Invalid parameterization"))

    ## Initialize objects to be passed to the output
    post2draws <- NA
    phi_pop <- NA
    stage2_summary <- NA
    
    if (multistage == TRUE) {
        ####################################################
        ## Stage 1: Run mlVAR model to extract residuals. ##
        ####################################################
        
        ## Model returns a list containing the phi pop vlaues random efects.
        ## Residuals are resid
        cat("\n---------------------------\n" )
        cat("Stage 1: mlVAR is estimated \n")
        cat("---------------------------\n" )
        fit_stage1 <- stage1_model$variational(data = stan_data,
                                               iter = 50000,
                                               threads = threads_per_chain,
                                               ...)
        ## Extract relevant params and residuals

        T <- stan_data$T
        nt <- stan_data$nt

        ## Extract Residuals and average across draws
        draws <- fit_stage1$draws(format = "draws_df")

        ## Phi0 population intercept
        phi0_pop <- suppressWarnings(colMeans(draws[, grep("^phi0_pop\\[", names(draws))]))
        ## corresponding SD's
        phi0_pop_sd <- suppressWarnings(apply(draws[, grep("^phi0_pop\\[", colnames(draws))], 2, sd, na.rm = TRUE))
        ## RE estimate:
        sigma_phi0_log <-  suppressWarnings(colMeans(draws[, grep("^sigma_phi0_log\\[", names(draws))]))
        sigma_phi0_log_sd <-  suppressWarnings(apply(draws[, grep("^sigma_phi0_log\\[", colnames(draws))], 2, sd, na.rm = TRUE))
        
        
        ## Pop phi:
        vec_phi_pop <- suppressWarnings(colMeans(draws[, grep("^vec_phi_pop\\[", names(draws))]))
        phi_pop <- matrix(vec_phi_pop, nt, nt)
        
        ## Subject-specific VAR matrices (posterior mean)
        extract_matrix <- function(prefix, J, nt) {
            arr <- array(NA_real_, dim = c(nt, nt, J))
            for (j in 1:J) {
                for (r in 1:nt) {
                    for (c in 1:nt) {
                        name <- sprintf("%s[%d,%d,%d]", prefix, j, r, c)  # e.g., "Phi_j[1,2,3]"
                        if (name %in% colnames(draws)) { arr[r, c, j] <- mean(draws[[name]])
                        }
                    }
                }
            }
            arr
        }

        ## Individual phi matrices: Not needed - commented out for now
        ## Phi_j_hat <- extract_matrix("Phi", J = J, nt = nt)

        ## phi0_j_hat <- array(NA_real_, dim = c(J, nt))
        ## for (j in 1:J) {
        ##     phi0_j_hat[j, ] <-  colMeans(draws[, grep(sprintf("^phi0_j\\[%d,", j), names(draws))])
        ## }

        ## Sigma (posterior mean)
        Sigma_hat <- matrix(NA_real_, nrow = nt, ncol = nt)
        for (r in 1:nt) for (c in 1:nt) {
                            nm <- sprintf("Sigma[%d,%d]", r, c)
                            if (nm %in% colnames(draws)) Sigma_hat[r, c] <- mean(draws[[nm]])
                        }

        ## Residuals: list J of T x nt matrices
        resid_list <- vector("list", J)
        for (j in 1:J) {
            mat <- matrix(NA_real_, nrow = T, ncol = nt)
            for (t in 1:T) {
                for (d in 1:nt) {
                    nm <- sprintf("resid[%d,%d,%d]", j, d, t)
                    if (nm %in% colnames(draws)) mat[t, d] <- mean(draws[[nm]])
                }
            }
            resid_list[[j]] <- mat
        }

        ## collect residuals and add back variable names
        residuals <- list()
        for (i in 1:J) {
            residuals[[i]] <- resid_list[[i]]
            colnames(residuals[[i]]) <- colnames(stan_data$rts[[1]])
        }

        ## Save residuals to pass along as data for other models
        stan_data$residuals <- residuals

        ## Add phi and Sigma to stan_data
        stan_data$phi0_pop_stage1 <- phi0_pop
        stan_data$phi0_pop_stage1_sd <- phi0_pop_sd
        stan_data$sigma_phi0_log_stage1 <- sigma_phi0_log
        stan_data$sigma_phi0_log_stage1_sd <- sigma_phi0_log_sd

        stan_data$phi_pop <- phi_pop

        ## stage1_params <- list(
        ##     phi0_pop = phi0_pop,
        ##     sigma_phi0_log = sigma_phi0_log
        ## )
        
        ###############################################################
        ## Stage 2:                                                  ##
        ##   Step 1: Compute random effects for uncond corr matrix S ##
        ##   Step 2: Compute hierarchical GARCH on diag of DRD       ##
        ###############################################################

        ## Step 1:
        ## Compute individual sample correlations
        ## Extract fisher z transform off diagional elements
        zhat <- lapply(seq_len(J), function(x) {
            Rj <- cor(residuals[[x]])
            atanh(Rj[lower.tri(Rj)])
        })
        ## S data
        Sdim <-  nt * (nt-1) / 2
        s_data <- list(
            J = J,
            nt = stan_data$nt,
            Sdim = Sdim,
            zhat = zhat
        )
        cat("\n--------------------------------------------\n" )
        cat("Stage 2, Step 1: Estimate random effect of S\n")
        cat("--------------------------------------------\n\n" )
        precomp_fit <- stage2_model$variational(data = s_data,
                                                iter = 50000,
                                                threads = threads_per_chain,
                                                ...)
        
        ## Extract the random effects SD vector
        post2draws <- precomp_fit$draws(format = "matrix", variables = paste0("sigma_z[", 1:Sdim, "]"))
        ## Simplest approach: Average across draws and use E(sigma_z) as SD for all S ranefs
        sigma_z_hat <- colMeans(post2draws)
        ## Attach to stan_data
        stan_data$Sdim <- Sdim
        stan_data$S_vec_tau_fixed <- sigma_z_hat

        ## Step 2:
        ## Use residuals from Stage 1
        cat("\n------------------------------------------\n" )
        cat("Stage 2, Step 2: Estimate GARCH on diag(D)\n")
        cat("------------------------------------------\n\n" )

        garch_data <- list(J = J,
                           nt = stan_data$nt,
                           T = T,
                           u = residuals)
        garch_fit <- stage2garch_model$variational(data = garch_data,
                                                   iter = 50000,
                                                   threads = threads_per_chain,
                                                   ...)
        ## Extract paramaters for downstream mlVAR-DCC model
        #### Garch h params
        draws_df <- garch_fit$draws(format = "draws_df")

        ## population fixed effects (directly use the raw / logit-scale)
        ## mu_c is on log scale
        mu_c <-     suppressWarnings(colMeans(draws_df[, grep("^mu_c\\[", colnames(draws_df))]))
        mu_a_raw <- suppressWarnings(colMeans(draws_df[, grep("^mu_a_raw\\[", colnames(draws_df))]))
        mu_b_raw <- suppressWarnings(colMeans(draws_df[, grep("^mu_b_raw\\[", colnames(draws_df))]))
        ## corresponding sampling SD's
        sd_mu_c     <- suppressWarnings(apply(draws_df[, grep("^mu_c\\[",     colnames(draws_df))], 2, sd, na.rm = TRUE))
        sd_mu_a_raw <- suppressWarnings(apply(draws_df[, grep("^mu_a_raw\\[", colnames(draws_df))], 2, sd, na.rm = TRUE))
        sd_mu_b_raw <- suppressWarnings(apply(draws_df[, grep("^mu_b_raw\\[", colnames(draws_df))], 2, sd, na.rm = TRUE))

        ## random effect SDs
        tau_c_log <- suppressWarnings(colMeans(draws_df[, grep("^c_h_tau_log\\[", colnames(draws_df))]))
        tau_a_log <- suppressWarnings(colMeans(draws_df[, grep("^a_h_tau_log\\[", colnames(draws_df))]))
        tau_b_log <- suppressWarnings(colMeans(draws_df[, grep("^b_h_tau_log\\[", colnames(draws_df))]))
        ## Correspondings sampling SD's
        sd_tau_c_log <- suppressWarnings(apply(draws_df[, grep("^c_h_tau_log\\[", colnames(draws_df))], 2, sd, na.rm = TRUE))
        sd_tau_a_log <- suppressWarnings(apply(draws_df[, grep("^a_h_tau_log\\[", colnames(draws_df))], 2, sd, na.rm = TRUE))
        sd_tau_b_log <- suppressWarnings(apply(draws_df[, grep("^b_h_tau_log\\[", colnames(draws_df))], 2, sd, na.rm = TRUE))


        ## Add to stan_data, all nt length vectors to pass to stan
        stan_data$c_h_fixed_s2 <- mu_c      # log-scale
        stan_data$a_h_fixed_s2 <- mu_a_raw  # logit-scale
        stan_data$b_h_fixed_s2 <- mu_b_raw  # logit-scale
        stan_data$c_h_fixed_s2_sd <- sd_mu_c      # log-scale
        stan_data$a_h_fixed_s2_sd <- sd_mu_a_raw  # logit-scale
        stan_data$b_h_fixed_s2_sd <- sd_mu_b_raw  # logit-scale

        stan_data$c_h_tau_log_s2 <- tau_c_log
        stan_data$a_h_tau_log_s2 <- tau_a_log
        stan_data$b_h_tau_log_s2 <- tau_b_log
        stan_data$c_h_tau_log_s2_sd <- sd_tau_c_log
        stan_data$a_h_tau_log_s2_sd <- sd_tau_a_log
        stan_data$b_h_tau_log_s2_sd <- sd_tau_b_log

        ## Assemble for easy read-out from fitted object
        stage2_summary <- list(
            c_h_fixed = mu_c,     # log-scale
            a_h_fixed = mu_a_raw, # logit-scale
            b_h_fixed = mu_b_raw, # logit-scale
            c_h_tau   = exp(tau_c_log),
            a_h_tau   = exp(tau_a_log),
            b_h_tau   = exp(tau_b_log)
        )

        cat("\n---------------------------------------\n")
        cat("\n Stage 3: Estimate final mlVAR-DCC\n")
        cat("\n---------------------------------------\n")
    }

    ## HMC Sampling
    ## If meanstructure = constant, only GARCH model is estimated on residuals from stge 1,
    ## otherwise data remains original TS and meanstructure is estimated as well
    if(meanstructure == 'constant') {
        stan_data$rts <- residuals
    }


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

    ## For mutlistage: Paste the post2draws from stage 2 to the fitted model 
    ## model_fit <- cbind(model_fit, post2draws)
    
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
                     S_pred = S_pred,
                     S_vec_tau_post = post2draws,
                     phi_pop = phi_pop,
                     stage2_summary = stage2_summary)

  class(return_fit) <- "dcnet"
  return(return_fit)
}


#' Models supported by dcnet
#'
#' To be used when checking whether a parameterization or object type is a supported type.
#' May facilitate more parameterizations, as we only have to update these, and the switch statements.
#' @keywords internal
#' @author Philippe Rast
supported_models <- c("CCC", "DCC", "DCCj", "DCCms")
