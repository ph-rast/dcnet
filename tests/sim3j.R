# ─────────────────────────────────────────────────────────────────────────────
# Full R script: future + progressr + possibly for 100 replications
# ─────────────────────────────────────────────────────────────────────────────

# 0) Install needed packages if not already installed

## Only install after repo update
#if (!requireNamespace("remotes")) { 
#  install.packages("remotes")   
#}   
#remotes::install_github("ph-rast/dcnet")

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

safe_sample <- function(s, replication_data, max_retries = 3) {
  fit_init <- dcnet(
    data             = replication_data[[1]],
    parameterization = "DCCj",
    J                = replication_data$N,
    group            = replication_data[[2]],
    standardize_data = FALSE,
    init             = 0,
    meanstructure    = "VAR",
    iterations       = 50000,
    tol_rel_obj      = 0.007,
    sampling_algorithm = "variational"
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
      "All proposed step-sizes failed|algorithm may not have converged|Exception:",
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
    'a_q_pop', 'l_a_q_sigma', 'b_q_pop', 'l_b_q_sigma', 'S_vec_fixed', 'S_vec_tau'
)
var <- c(
    "phi0_fixed", "phi0_sd", "fixed_phi", "phi_ranef_sd", "phi_ranef_sd",
    "phi_sd_diag", "phi_sd_off",
    "log_c_fixed", "log_c_r_sd", "a_h_fixed", "a_h_r_sd", "b_h_fixed", "b_h_r_sd",
    "l_a_q_fixed", "l_a_q_r_sd", "l_b_q_fixed", "l_b_q_r_sd", "fixed_S_atanh", "ranS_sd"
)
stopifnot(length(variables_m) == length(var))

# 5) Build condition grid + define 100 replications
ns        <- c(25, 50, 100)
tls       <- c(25, 50, 100)
simcond   <- expand.grid(N = ns, tl = tls)
n_conds   <- nrow(simcond)
n_reps    <- 10
task_grid <- expand.grid(idx = seq_len(n_conds), rep = seq_len(n_reps))

# 6) Setup future + progressr
plan(multisession, workers = 56)
handlers("txtprogressbar")

# 7) Wrap core task in a “possibly” to swallow errors:
run_one <- possibly(function(idx, rep) {
  p()  # tick the progress bar
  N  <- simcond$N[idx]
  tl <- simcond$tl[idx]

  dat <- simulate_data(N = N, tl = tl, nts = 3)
  fit <- safe_sample(rep, replication_data = dat)
  if (is.null(fit)) return(NULL)

  out <- lapply(seq_along(variables_m), function(pi) {
    stan_par <- variables_m[pi]
    sim_var  <- var[pi]

    draws_mat <- fit$model_fit$draws(
      variables = stan_par, format = "matrix"
    )
    truth_vec <- unlist(dat[[3]][[sim_var]])

    SF <-
      ##if (stan_par == "S_vec_tau") {
      ##tmp <- t(apply(fit$S_vec_tau_post, 2,
      ##               quantile, probs = c(0.5,0.025,0.975)))
      ##colnames(tmp) <- c("mean","q2.5","q97.5")
      ##as.data.frame(tmp)
    ##} else {
      fit$model_fit$summary(
        variables       = stan_par,
        fun             = "mean",
        extra_quantiles = ~ posterior::quantile2(., probs = c(.025, .975))
      )
    ##}

    data.frame(
      N           = N,
      tl          = tl,
      replication = rep,
      parameter   = stan_par,
      coverage    = mean(truth_vec > SF$q2.5 & truth_vec < SF$q97.5),
      bias        = mean(as.numeric(draws_mat) -
                           rep(truth_vec, each = nrow(draws_mat))),
      rmse        = sqrt(mean((as.numeric(draws_mat) -
                                 rep(truth_vec, each = nrow(draws_mat)))^2)),
      ci_width    = mean(abs(SF$q97.5 - SF$q2.5)),
      looic       = looic(fit),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}, otherwise = NULL)

# 8) Run with live progress
with_progress({
  p <- progressor(along = seq_len(nrow(task_grid)))
  per_rep <- future_pmap_dfr(
    .l       = task_grid,
    .f       = run_one,
    .options = furrr_options(seed = TRUE)
  )
})

# 9) Aggregate over 100 reps
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
#saveRDS(final_results, "simulation_performance_100reps.rds")

final_results <- readRDS("simulation_performance_j_100reps.rds")

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
