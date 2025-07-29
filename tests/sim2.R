######################################
## Simulate 
######################################
devtools::load_all()

## Simulation conditions: Only vary sample size and time series length
ns <- c(30, 60, 100)
tls <- c(30, 60, 100)

simcond <- expand.grid(N = ns, tl = tls)
simcond



## list containing all replicaion datasets:
replication_data <- list()

## Create data:
simulate_data <- function(N = 10, tl = 5, nts = 3) {
    simdat <- dcnet:::.simuDCC(
        tslength = tl, N = N, n_ts = nts,
        phi_mu = 0, ## populate phi
        phi0_sd = 0.4, ## create variation in the intercepts of the time series
        phi_sd_diag = 0.3,  ## SD of own lags across TS's
        phi_sd_off = 0.05,  ## SD of cross lags across TS's
        phi_ranef_sd = 0.01, ## random effects in phi
        log_c_fixed = rep(-0.9, nts),
        log_c_r_sd = 0.3,
        a_h_fixed = rep(-2.5, nts),
        a_h_r_sd = 0.1,
        b_h_fixed = rep(-1.5, nts),  ## On logit scale
        b_h_r_sd = 0.1,
        l_a_q_fixed = -1.5,  ## on logit scale
        l_b_q_fixed = -.5,   ## on logit scale
        l_a_q_r_sd = 0.2,
        l_b_q_r_sd = 0.2,
        phi0_fixed =  rep(0, nts),
        ranS_sd = 0.25, ## random effects on atanh scale
        stationarity_phi = FALSE
    )
    rtsgen <- lapply(seq(dim(simdat[[1]])[3]), function(x) t(simdat[[1]][, , x]))
    groupvec <- rep(c(1:N), each = tl)
    return(list(rtsgen, groupvec, simdat, N = N))
}


## Function to prevent loop to stop when stan model encounters errors
## Try fitting model 3 times before moving on

safe_sample <- function(s, max_retries = 3, replication_data) {
    replication_data <- replication_data

    ## First run meanfield algo to obtain init values:
    fit_init <- dcnet(
        data = replication_data[[1]], parameterization = "DCCrs", J = replication_data$N,
        group = replication_data[[2]], standardize_data = FALSE,
        init = 0,
        meanstructure = "VAR",
        iterations = 50000,
        # eta = 0.05,
        tol_rel_obj = 0.005,
        sampling_algorithm = "variational"
    )
    
    ## Helper function to check if model_fit is broken
    is_model_fit_broken <- function(fit) {
        tryCatch({
            if (!inherits(fit$model_fit, "CmdStanFit")) return(TRUE)
            fit$model_fit$metadata()  # harmless method, fails if object is invalid
            return(FALSE)
        }, error = function(e) {
            return(TRUE)
        })
    }

    ## Fit model with fullrank
    for (attempt in seq_len(max_retries)) {
        message(sprintf("Fitting attempt %d for index %d", attempt, s))

        fit <- tryCatch({
            fit_init
        }, error = function(e) {
            message(sprintf("Hard failure on attempt %d: %s", attempt, e$message))
            return(NULL)
        })

        if (is.null(fit)) next

        ## Check if model_fit is broken
        if (is_model_fit_broken(fit)) {
            message("model_fit is invalid or unusable — retrying...")
            next
        }

        ## Try to extract model output
        output_lines <- tryCatch({
            out <- fit$model_fit$output()
            if (is.character(out)) out else as.character(unlist(out))
        }, error = function(e) {
            message("Could not retrieve model output — assuming failure.")
            return("ERROR_FETCHING_OUTPUT")
        })

        ## Patterns that indicate failure or non-convergence
        failure_patterns <- c(
            "All proposed step-sizes failed",
            "algorithm may not have converged",
            "Exception:.*not positive definite",
            "Exception:",
            "Fitting failed. Unable to print"
        )

        failed <- any(sapply(failure_patterns, function(pat) {
            any(grepl(pat, output_lines, ignore.case = TRUE, perl = TRUE))
        }))

        if (failed) {
            message("Detected Stan warning or internal failure — retrying...")
            next
        }

        ## Everything looks okay
        return(fit)
    }
    message(sprintf("All %d attempts failed for index %d.", max_retries, s))
    return(NULL)
}


## Metrics: 
## Function to compute coverage probability ranging from L to U in the original population distribution
overlap <- function(population, L, U) {
    coverage <- as.numeric(population > L & population < U)
    return(coverage)
}

crirange <- function(L, U) {
    bwidth <- abs(as.numeric(U-L))
    return(bwidth)
}


rmse <- function(model, population) {
    rmse <- sqrt(mean((model - population)^2))
    return(rmse)
}

bias <- function(model, population) {
    bias <- mean(model - population)
    return(bias)
}

sbc <- function(model, population, column) {
  ## Algorithm 1 from: Validating Bayesian Inference Algorithms with Simulation-Based Calibration (2020)
  bin <- list()
  binl <- 0
  ## Assuming 1000 draws in the fit objects (standard)
  ## This creates 20 bins 
  for(start in seq(1, 951, by=50)) {
    binl <- binl+1
    end <- start + 49
    current_sequence <- start:end
    bin[binl] <- sum(model[start:end,column] < population[start:end,column])
  }
  return(bin)
}

looic <- function(fit) {
    log_lik_r <- grep( "log_lik", colnames(fit$model_fit$draws( )) )
    log_lik_r <- fit$model_fit$draws( )[, log_lik_r]
    r_eff_r <- loo::relative_eff(exp(log_lik_r ),  chain_id = rep(1,  1000))
    fr <- loo::loo(log_lik_r, r_eff = r_eff_r)
    looic <- fr$estimates[3,1]
    return(looic)
}

## ## Stan variables:
## variables_m <- c(
##     'phi0_fixed', 'phi0_tau', 'vec_phi_fixed', 'sigma_re_own', 'sigma_re_cross',
##     'tau_own_log', 'tau_cross_log',
##     'c_h_fixed', 'c_h_tau', 'a_h_fixed', 'a_h_tau', 'b_h_fixed', 'b_h_tau',
##     'l_a_q', 'l_a_q_sigma', 'l_b_q', 'l_b_q_sigma', 'S_vec_fixed', 'S_vec_tau')

## ## Simulation variables: Note that the simulatin script only defines one random
## ## effect for phi, phi_ranef_sd, but the stan model captures the random effects
## ## in both, the own and cross lag of phi in sigma_re_own/cross
## var <- c(
##     "phi0_fixed", "phi0_sd", "fixed_phi", "phi_ranef_sd", "phi_ranef_sd",
##     "phi_sd_diag", "phi_sd_off",
##     "log_c_fixed", "log_c_r_sd", "a_h_fixed", "a_h_r_sd", "b_h_fixed", "b_h_r_sd",
##     "l_a_q_fixed", "l_a_q_r_sd", "l_b_q_fixed", "l_b_q_r_sd", "fixed_S_atanh", "ranS_sd"
## )

## if(!length(variables_m) == length(var)) stop("variable list does not match!")

## # Initialize empty lists to store the results for each variable
## cov_list <- list()
## crir_list <- list()
## rmse_list <- list()
## bias_list <- list()
## bins_list <- list()
## looic_list <- list()

## for (s in 1:10) {

##     replication_data <- simulate_data(N = simcond$N[1], tl = simcond$tl[1])
##     fit_r <- safe_sample(s, replication_data = replication_data)

##     if (is.null(fit_r)) {
##         next
##     }

##     for (p in 1:length(variables_m)) {
##         SF <- fit_r$model_fit$summary(
##             variables = variables_m[p], "mean",
##             extra_quantiles = ~ posterior::quantile2(., probs = c(0.025, 0.975))
##         )
##         cov_list[[variables_m[p]]][[s]] <-
##             overlap(
##                 unlist(replication_data[[3]][var[p]]),
##                 L = SF$q2.5, U = SF$q97.5
##             )
##         crir_list[[variables_m[p]]][[s]] <-
##             crirange(L = SF$q2.5, U = SF$q97.5)
##         rmse_list[[variables_m[p]]][[s]] <-
##             sapply(seq_len(nrow(SF)), function(i) {
##                 rmse(
##                     fit_r$model_fit$draws(variables = variables_m[p])[, i],
##                     unlist(replication_data[[3]][var[p]])
##                 )
##             })
##         bias_list[[variables_m[p]]][[s]] <-
##             sapply(seq_len(nrow(SF)), function(i) {
##                 bias(
##                     fit_r$model_fit$draws(variables = variables_m[p])[, i],
##                     unlist(replication_data[[3]][var[p]])
##                 )
##             })
##         looic_list[[variables_m[p]]][[s]] <-
##             suppressWarnings(looic(fit_r))
##     }
##     gc()
##     s
## }

## ## check raw lists:
## cov_list
## crir_list
## rmse_list
## bias_list
## bins_list
## looic_list

## ## compute averages across list
## var_averages <- list()
## for (i in 1:length(cov_list)) {
##     nested_average <- c()                                             
##     for (j in 1:length(cov_list[[i]])) {
##         if (is.null(cov_list[[i]][[j]])) cov_list[[i]][[j]] <- NA
##         nested_average <- c(nested_average, mean(cov_list[[i]][[j]], na.rm = TRUE))
##         var_averages[[i]] <- nested_average
##     }
## }
## names(var_averages) <- names(cov_list)
## var_averages


## var_av <- sapply(var_averages, function(x ) mean(x, na.rm = TRUE))
## var_av


## bias_averages <- list()
## for (i in 1:length(bias_list)) {
##     nested_average <- c()
##     for (j in 1:length(bias_list[[i]])) {
##         if (is.null(bias_list[[i]][[j]])) bias_list[[i]][[j]] <- NA
##         nested_average <- c(nested_average, mean(bias_list[[i]][[j]], na.rm = TRUE))
##         bias_averages[[i]] <- nested_average
##     }
## }
## names(bias_averages) <- names(bias_list)
## bias_averages

## bias_av <- sapply(bias_averages, function(x) mean(x, na.rm = TRUE))
## round(cbind(bias_av, var_av), 3)



## n30tl50 <- round(cbind(bias_av, var_av), 3)
## n30tl50
## ##n50tl50 <- round(var_av, 2)
## n50tl50
## ## n50tl75 <- round(var_av, 2)
## n50tl75
## ## n75tl75 <- round(var_av, 2)
## n75tl75
## ## n75tl100 <- round(var_av, 2)
## n75tl100
## n100tl100 <- round(cbind(bias_av, var_av), 3)
## n100tl100



### PArallel version:
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   
remotes::install_github("ph-rast/dcnet")


# at the top of your script:
library(doParallel)
library(foreach)
library(dcnet) ## needs to be loaded for parallel computation
#library(loo)       # for looic()
#library(posterior) # for as_draws_matrix()

# set up your 9 conditions:
ns  <- c(25, 50, 100)
tls <- c(25, 50, 100)
simcond <- expand.grid(N = ns, tl = tls)
n_conds <- nrow(simcond)
replications <- 20

# register one worker per condition:
cl <- makeCluster(min(n_conds, 28))
registerDoParallel(cl)

# define the Stan ↔ sim parameter names (positional matching):
## Stan variables:
variables_m <- c(
    'phi0_fixed', 'phi0_tau', 'vec_phi_fixed', 'sigma_re_own', 'sigma_re_cross',
    'tau_own_log', 'tau_cross_log',
    'c_h_fixed', 'c_h_tau', 'a_h_fixed', 'a_h_tau', 'b_h_fixed', 'b_h_tau',
    'l_a_q', 'l_a_q_sigma', 'l_b_q', 'l_b_q_sigma', 'S_vec_fixed', 'S_vec_tau')

## Simulation variables: Note that the simulatin script only defines one random
## effect for phi, phi_ranef_sd, but the stan model captures the random effects
## in both, the own and cross lag of phi in sigma_re_own/cross
var <- c(
    "phi0_fixed", "phi0_sd", "fixed_phi", "phi_ranef_sd", "phi_ranef_sd",
    "phi_sd_diag", "phi_sd_off",
    "log_c_fixed", "log_c_r_sd", "a_h_fixed", "a_h_r_sd", "b_h_fixed", "b_h_r_sd",
    "l_a_q_fixed", "l_a_q_r_sd", "l_b_q_fixed", "l_b_q_r_sd", "fixed_S_atanh", "ranS_sd"
)

stopifnot(length(variables_m)==length(var))

# Make the local functions & objects available on each worker
clusterExport(cl, c("simulate_data", "safe_sample", "variables_m", "var", "simcond"))

# Load packages in each worker
clusterEvalQ(cl, {
  library(dcnet)
  library(loo)
  library(posterior)
  TRUE
})

# run each condition in parallel:
results_list <- foreach(i = seq_len(n_conds),
                        .packages = c("dcnet","loo","posterior")) %dopar% {
  N  <- simcond$N[i]
  tl <- simcond$tl[i]

  #--- replicate the condition:
  cov_mat  <- matrix(NA, nrow = length(variables_m), ncol = replications)
  bias_mat <- matrix(NA, nrow = length(variables_m), ncol = replications)
  rmse_mat <- matrix(NA, nrow = length(variables_m), ncol = replications)
  crir_mat <- matrix(NA, nrow = length(variables_m), ncol = replications)
  loo_vec  <- numeric(replications)

  for (s in seq_len(replications)) {
    # 1) simulate
    rep_data <- simulate_data(N = N, tl = tl, nts = 3)
    # 2) fit (meanfield VI)
    fit_r <- safe_sample(s, replication_data = rep_data)
    if (is.null(fit_r)) next

    # 3) extract draws & summaries
    for (p in seq_along(variables_m)) {
      stan_par <- variables_m[p]
      sim_name <- var[p]

      # posterior draws as matrix [iter x length(truth)]
      draws_mat <- fit_r$model_fit$draws(
                                       variables = stan_par,
                                       format = "matrix")
      truth_vec <- unlist(rep_data[[3]][[sim_name]])

      # summary with 2.5 & 97.5 % quantiles
      SF <- fit_r$model_fit$summary(
        stan_par, "mean",
        extra_quantiles = ~ posterior::quantile2(., probs = c(0.025, 0.975))
      )

      # coverage for each element, then average
      cov_mat[p, s]  <- mean(truth_vec > SF$q2.5 & truth_vec < SF$q97.5)

      # bias and RMSE (flatten draws & truth)
      all_draws <- as.numeric(draws_mat)
      all_truth <- rep(truth_vec, each = nrow(draws_mat))
      bias_mat[p, s] <- mean(all_draws - all_truth)
      rmse_mat[p, s] <- sqrt(mean((all_draws - all_truth)^2))

      # average 95% CI width (elementwise)
      crir_mat[p, s] <- mean(abs(SF$q97.5 - SF$q2.5))
    }

    # 4) LOOIC for this replicate
      loo_vec[s] <- looic(fit_r)
  }

  # 5) Aggregate across 10 replications
  data.frame(
    N        = N,
    tl       = tl,
    parameter= variables_m,
    coverage = rowMeans(cov_mat, na.rm = TRUE),
    bias     = rowMeans(bias_mat, na.rm = TRUE),
    rmse     = rowMeans(rmse_mat, na.rm = TRUE),
    crir     = rowMeans(crir_mat, na.rm = TRUE),
    looic    = mean(loo_vec, na.rm = TRUE)
  )
}

stopCluster(cl)

# bind into one big data.frame:
final_results <- do.call(rbind, results_list)

head(final_results)
# save for later:
saveRDS(final_results, "simulation_performance.rds")


## Plots:

##reshape for plotting
library(tidyverse)

fr_long <- final_results %>%
  pivot_longer(
    cols      = coverage:looic,
    names_to  = "metric",
    values_to = "value"
  )

fr_long

## coverage heatmap
fr_long %>%
  filter(metric == "coverage") %>%
  ggplot(aes(x = factor(N), y = factor(tl), fill = value)) +
    geom_tile(color = "white") +
    facet_wrap(~ parameter, ncol = 3) +
    scale_fill_viridis_c(limits = c(0,1)) +
    labs(
      x = "N (subjects)",
      y = "T (timepoints)",
      fill = "Coverage",
      title = "95% CI Coverage by Condition and Parameter"
    ) +
    theme_minimal(base_size = 12)

## Bias
fr_long %>%
  filter(metric == "bias") %>%
  ggplot(aes(x = N, y = value, color = factor(tl))) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
    labs(
      x = "N (subjects)",
      y = "Bias",
      color = "T (timepoints)"
    ) +
  theme_minimal()

# 2. Create a named list of ggplot objects
plot_list <- list(
  coverage = fr_long %>%
    filter(metric == "coverage") %>%
    ggplot(aes(factor(N), factor(tl), fill = value)) +
      geom_tile(color = "white") +
      facet_wrap(~ parameter, ncol = 3) +
      scale_fill_viridis_c(limits = c(0,1)) +
      labs(x = "N", y = "T", fill = "Coverage") +
      theme_minimal(),

  bias = fr_long %>%
    filter(metric == "bias") %>%
    ggplot(aes(x = N, y = value, color = factor(tl))) +
      geom_line() + geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
      labs(x = "N", y = "Bias", color = "T") +
      theme_minimal(),

  rmse = fr_long %>%
    filter(metric == "rmse") %>%
    ggplot(aes(x = N, y = value, color = factor(tl))) +
      geom_line() +
      facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
      labs(x = "N", y = "RMSE", color = "T") +
      theme_minimal(),

  ci_width = fr_long %>%
    filter(metric == "crir") %>%
    ggplot(aes(x = N, y = value, color = factor(tl))) +
      geom_line() +
      facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
      labs(x = "N", y = "CI width", color = "T") +
      theme_minimal(),

  looic = fr_long %>%
    filter(metric == "looic") %>%
    ggplot(aes(x = N, y = value, color = factor(tl))) +
      geom_line() +
      facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
      labs(x = "N", y = "LOOIC", color = "T") +
      theme_minimal()
)

# 3. Write them all to a multi‐page PDF
pdf("simulation_diagnostics.pdf", width = 11, height = 8.5)
for (nm in names(plot_list)) {
  print(plot_list[[nm]] + ggtitle(str_to_title(nm)))
}
dev.off()

# 4. (Optional) save the list of plots for later use
saveRDS(plot_list, "plot_list.rds")

## table
summary_table <- final_results %>%
  group_by(parameter) %>%
  summarise(
    cov_min   = min(coverage),
    cov_mean  = mean(coverage),
    bias_abs  = mean(abs(bias)),
    rmse_mean = mean(rmse),
    ciw_mean  = mean(crir),
    looic     = mean(looic)
  ) %>%
  arrange(cov_mean)
print(summary_table)
