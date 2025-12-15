# ─────────────────────────────────────────────────────────────────────────────
# Full R script: future + progressr + possibly for 100 replications
# ─────────────────────────────────────────────────────────────────────────────

# 0) Install needed packages if not already installed

## Only install after repo update
#if (!requireNamespace("remotes")) { 
#  install.packages("remotes")   
#}   
#remotes::install_githlub("ph-rast/dcnet")

needed <- c("furrr", "progressr", "purrr", "loo", "posterior",
            "dplyr", "tidyr")
to_install <- needed[!needed %in% installed.packages()[, "Package"]]
if (length(to_install)) install.packages(to_install)

# 1) Load libraries
library(furrr)
library(progressr)
library(purrr)
library(dcnet)      # for simulate_data() & safe_sample()
library(loo)        # for loo & relative_eff
library(posterior)
library(dplyr)
library(tidyr)

# 2) Define looic() helper
looic <- function(fit) {
    draws_mat <- fit$model_fit$draws(format = "matrix")
    ll_cols   <- grep("^log_lik", colnames(draws_mat))
    log_lik   <- draws_mat[, ll_cols, drop = FALSE]
    r_eff     <- loo::relative_eff(exp(log_lik),
                                   chain_id = rep(1, nrow(log_lik)))
    loo_res   <- loo::loo(log_lik, r_eff = r_eff)
    loo_res$estimates["looic", "Estimate"]
}

# 3) Your helper functions (copy in your actual versions)
simulate_data <- function(N = 10, tl = 5, nts = 3) {
  simdat   <- dcnet:::.simuDCC(
    tslength         = tl, N = N, n_ts = nts,
    phi_mu           = 0, phi0_sd     = 0.4,
    phi_sd_diag      = 0.3, phi_sd_off = 0.05,
    phi_ranef_sd     = 0.01,
    log_c_fixed      = rep(-0.9, nts), log_c_r_sd = 0.3,
    a_h_fixed        = rep(-2.5, nts), a_h_r_sd  = 0.1,
    b_h_fixed        = rep(-1.5, nts), b_h_r_sd  = 0.1,
    l_a_q_fixed      = -1.5, l_b_q_fixed = -0.5,
    l_a_q_r_sd       = 0.2, l_b_q_r_sd = 0.2,
    phi0_fixed       = rep(0, nts), ranS_sd    = 0.25,
    stationarity_phi = FALSE
  )
  rtsgen   <- lapply(seq(dim(simdat[[1]])[3]),
                     function(x) t(simdat[[1]][,,x]))
  groupvec <- rep(seq_len(N), each = tl)
  list(rtsgen, groupvec, simdat, N = N)
}

safe_sample <- function(s, replication_data, max_retries = 2) {
  fit_init <- dcnet(
    data             = replication_data[[1]],
    parameterization = "DCCms",
    J                = replication_data$N,
    group            = replication_data[[2]],
    standardize_data = FALSE,
    init             = 0,
    meanstructure    = "VAR",
    iterations       = 1500,
    chains = 1,
    sampling_algorithm = "hmc"
  )

  is_broken <- function(f) {
    tryCatch({
      if (!inherits(f$model_fit, "CmdStanFit")) return(TRUE)
      f$model_fit$metadata(); FALSE
    }, error = function(e) TRUE)
  }

  for (at in seq_len(max_retries)) {
    if (is_broken(fit_init)) next
    out_lines <- tryCatch({
      o <- fit_init$model_fit$output()
      as.character(unlist(o))
    }, error = function(e) "ERROR")
    fails <- any(grepl(
      "ERROR|Error|finished unexpectedly",
      out_lines, ignore.case = TRUE
    ))
    if (!fails) return(fit_init)
  }
  NULL
}

# 4) Stan ↔ sim variable mappings
variables_m <- c(
    'phi0_fixed', 'phi0_tau', 'vec_phi_fixed', 'sigma_re_own', 'sigma_re_cross',
    'tau_own', 'tau_cross',
    'c_h_fixed', 'c_h_tau', 'a_h_fixed', 'a_h_tau', 'b_h_fixed', 'b_h_tau',
    'l_a_q', 'l_a_q_sigma', 'l_b_q', 'l_b_q_sigma', 'S_vec_fixed', 'S_vec_tau'
)
var <- c(
    "phi0_fixed", "phi0_sd", "fixed_phi", "phi_ranef_sd", "phi_ranef_sd",
    "phi_sd_diag", "phi_sd_off",
    "log_c_fixed", "log_c_r_sd", "a_h_fixed", "a_h_r_sd", "b_h_fixed", "b_h_r_sd",
    "l_a_q_fixed", "l_a_q_r_sd", "l_b_q_fixed", "l_b_q_r_sd", "fixed_S_atanh", "ranS_sd"
)
stopifnot(length(variables_m) == length(var))

# 5) Build condition grid + define 100 replications
ns        <- c(50)#c(25,  50, 100)#, 150)
tls       <- c(50)#c(25, 50, 100)#, 150)
simcond   <- expand.grid(N = ns, tl = tls)
n_conds   <- nrow(simcond)
n_reps    <- 1
task_grid <- expand.grid(idx = seq_len(n_conds), reps = seq_len(n_reps))

## 6) Setup future + progressr
cores  <- parallelly::availableCores()
cores <- 4
chains <- 1 ## take chains from safe_sample function
workers <- max(1, floor((cores - 1) / chains))
# plan(multisession, workers = workers) ## Don't spawn here - do it later and remove workers after use (or feel the wrath of John... )
handlers("txtprogressbar")

# 7) Wrap core task in a “possibly” to swallow errors:
run_one <- possibly(function(idx, reps) {
    p()  # tick the progress bar
    N  <- simcond$N[idx]
    tl <- simcond$tl[idx]

    dat <- simulate_data(N = N, tl = tl, nts = 3)
    ## ## Recover per-person a_q, b_q used by the simulator to build Q_t
    ## ## a_q_j <- 1 / (1 + exp(-(dat[[3]]$l_a_q_fixed + rnorm(dat$N, 0, dat[[3]]$l_a_q_r_sd))))
    ## ## b_q_j <- (1 - a_q_j) * (1 / (1 + exp(-(dat[[3]]$l_b_q_fixed + rnorm(dat$N, 0, dat[[3]]$l_b_q_r_sd)))))
    ## a_q_j <- dat[[3]]$a_q
    ## b_q_j <- dat[[3]]$b_q
    ## ## Population means on probability scale:
    ## a_q_pop_true <- mean(a_q_j)
    ## b_q_pop_true <- mean(b_q_j)
    ## ## Truths on the *logit* scale to match Stan’s l_* parameters:
    ## l_a_q_truth <- qlogis(a_q_pop_true)
    ## l_b_q_truth <- qlogis(b_q_pop_true)
    ## ## fit <- safe_sample(reps, replication_data = dat)
    ## l_a_q_sigma_truth <- sd(qlogis(a_q_j))
    ## l_b_q_sigma_truth <- sd(qlogis(b_q_j))


    fit <- safe_sample(reps, replication_data = dat)
    if (is.null(fit)) {
        return(NULL)
    }

    out <- lapply(seq_along(variables_m), function(pi) {
        stan_par <- variables_m[pi]
        sim_var  <- var[pi]
        draws_mat <- fit$model_fit$draws(
                                       variables = stan_par, format = "matrix"
                                   )
        ## truth_vec <- switch(stan_par,
        ##                     "l_a_q"        = rep(l_a_q_truth,  length.out = ncol(draws_mat)),
        ##                     "l_b_q"        = rep(l_b_q_truth,  length.out = ncol(draws_mat)),
        ##                     ##"l_a_q_sigma"  = rep(dat[[3]]$l_a_q_r_sd, length.out = ncol(draws_mat)),
        ##                     ##"l_b_q_sigma"  = rep(dat[[3]]$l_b_q_r_sd, length.out = ncol(draws_mat)),
        ##                     "l_a_q_sigma"  = rep(l_a_q_sigma_truth, length.out = ncol(draws_mat)),
        ##                     "l_b_q_sigma"  = rep(l_b_q_sigma_truth, length.out = ncol(draws_mat)),
        ##                     ## default:
        ##                     unlist(dat[[3]][[sim_var]])
        ##                     )
        truth_vec <- unlist(dat[[3]][[sim_var]])
        SF <- if (stan_par == "S_vec_tau") {
                tmp <- t(apply(fit$S_vec_tau_post, 2,
                               quantile, probs = c(0.5,0.025,0.975)))
                colnames(tmp) <- c("mean","q2.5","q97.5")
                as.data.frame(tmp)
              } else {
                fit$model_fit$summary(
                                variables       = stan_par,
                                fun             = "mean",
                                extra_quantiles = ~ posterior::quantile2(., probs = c(.025, .975))
                              )
              }
        data.frame(
          N           = N,
          tl          = tl,
          replication = reps,
          parameter   = stan_par,
          coverage    = mean(truth_vec > SF$q2.5 & truth_vec < SF$q97.5),
          bias        = mean(as.numeric(draws_mat) -
                             rep(truth_vec, each = nrow(draws_mat))),
          rmse        = sqrt(mean((as.numeric(draws_mat) -
                                   rep(truth_vec, each = nrow(draws_mat)))^2)),
          ci_width    = mean(abs(SF$q97.5 - SF$q2.5)),
          looic       = suppressWarnings(looic(fit)),
          stringsAsFactors = FALSE
        )
    })

    do.call(rbind, out)
}, otherwise = NULL)

#  Run with live progress
## with_progress({
##     p <- progressor(steps = nrow(task_grid))
##     per_rep <- future_pmap_dfr(
##         .l       = task_grid,
##         .f       = run_one,
##         .options = furrr_options(seed = TRUE)
##     )
## })


##old_plan <-
future::plan(multisession, workers = workers)
##on.exit(future::plan(old_plan), add = TRUE)  # auto-restore when the block exits

with_progress({
  p <- progressr::progressor(steps = nrow(task_grid))
  per_rep <- furrr::future_pmap_dfr(
                      task_grid,
                      function(idx, reps) {
                        res <- run_one(idx, reps)  # do the heavy work first
                        p(sprintf("cond %d, rep %d done", idx, reps))  # tick on completion
                        res
                      },
                      .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)
                    )
})
future::plan(sequential)

per_rep

## 9) Aggregate over all reps
final_results <- per_rep %>%
  group_by(N, tl, parameter) %>%
  summarise(
    coverage = mean(coverage, na.rm = TRUE),
    bias     = mean(bias,     na.rm = TRUE),
    rmse     = mean(rmse,     na.rm = TRUE),
    ci_width = mean(ci_width, na.rm = TRUE),
    looic    = mean(looic,    na.rm = TRUE),
    .groups  = "drop"
  )

data.frame(final_results)

# 10) Save
# saveRDS(final_results, "simulation_performance_10reps.rds")

final_results <- readRDS("simulation_performance_ms_5reps.rds")
final_results <- readRDS("simulation_performance_20reps.rds")
final_results2 <- readRDS("simulation_performance_10Areps.rds")

head(final_results[,4:8])
final_result <- (final_results[, 4:8] + final_results2[,4:8])/2
final_result <- cbind(final_results[,1:3], final_result)

final_result

## Plots:

##reshape for plotting
library(tidyverse)

fr_long <- final_result %>%
  pivot_longer(
    cols      = coverage:looic,
    names_to  = "metric",
    values_to = "value"
  )

print(fr_long, n = 30)
dim(fr_long )

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
    filter(metric == "ci_width") %>%
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

## Write them all to a multi‐page PDF
pdf("hmc_simulation_diagnostics_10v2.pdf", width = 11, height = 8.5)
for (nm in names(plot_list)) {
  print(plot_list[[nm]] + ggtitle(str_to_title(nm)))
}
dev.off()

## save the list of plots for later use
saveRDS(plot_list, "hmc_plot_list.rds")

## table
summary_table <- final_results %>%
  group_by(parameter) %>%
  summarise(
    cov_min   = min(coverage),
    cov_mean  = mean(coverage),
    bias_abs  = mean(abs(bias)),
    rmse_mean = mean(rmse),
    ciw_mean  = mean(ci_width),
    looic     = mean(looic)
  ) %>%
  arrange(cov_mean)

print(summary_table)
