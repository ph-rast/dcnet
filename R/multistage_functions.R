##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param fit 
##' @param nt 
##' @param J 
##' @param T 
##' @return list
##' @keywords internal
##' @importFrom posterior as_draws_df
##' @author Philippe Rast
.extract_stage1_posterior <- function(fit, nt, J, T) {
  draws_df <- posterior::as_draws_df(fit$draws(format = "draws_df"))

  safe_mean <- function(pattern, expected_len = NULL) {
    cols <- grep(pattern, colnames(draws_df), perl = TRUE)
    if (length(cols) == 0) {
      warning("No columns matching pattern: ", pattern)
      if (!is.null(expected_len)) return(rep(NA_real_, expected_len))
      return(NA_real_)
    }
    out <- suppressWarnings(colMeans(draws_df[, cols, drop = FALSE]))
    if (!is.null(expected_len) && length(out) != expected_len) {
      warning("Length mismatch for ", pattern, ": got ", length(out), " expected ", expected_len)
    }
    out
  }

  # Population-level estimates
  phi0_pop_hat <- safe_mean("^phi0_pop\\[", expected_len = nt)
  tau_own_log_hat   <- if ("tau_own_log" %in% names(draws_df)) mean(draws_df$tau_own_log, na.rm = TRUE) else { warning("Missing tau_own_log"); NA }
  tau_cross_log_hat <- if ("tau_cross_log" %in% names(draws_df)) mean(draws_df$tau_cross_log, na.rm = TRUE) else { warning("Missing tau_cross_log"); NA }
  lambda_own_log_hat   <- safe_mean("^lambda_own_log\\[", expected_len = nt)
  lambda_cross_log_hat <- safe_mean("^lambda_cross_log\\[", expected_len = nt * (nt - 1))
  g_own_hat   <- safe_mean("^g_own\\[", expected_len = nt)
  g_cross_hat <- safe_mean("^g_cross\\[", expected_len = nt * (nt - 1))

  sigma_re_own_log_hat   <- if ("sigma_re_own_log" %in% names(draws_df)) mean(draws_df$sigma_re_own_log, na.rm = TRUE) else { warning("Missing sigma_re_own_log"); NA }
  sigma_re_cross_log_hat <- if ("sigma_re_cross_log" %in% names(draws_df)) mean(draws_df$sigma_re_cross_log, na.rm = TRUE) else { warning("Missing sigma_re_cross_log"); NA }
  sigma_phi0_hat <- if ("sigma_phi0" %in% names(draws_df)) mean(draws_df$sigma_phi0, na.rm = TRUE) else { warning("Missing sigma_phi0"); NA }

  # Reconstruct or grab Sigma
  Sigma_hat <- NULL
  if (any(grepl("^Sigma\\[1,1\\]", colnames(draws_df)))) {
    Sigma_hat <- matrix(NA_real_, nt, nt)
    for (r in 1:nt) for (c in 1:nt) {
      nm <- sprintf("Sigma[%d,%d]", r, c)
      if (nm %in% colnames(draws_df)) {
        Sigma_hat[r, c] <- mean(draws_df[[nm]])
      }
    }
  } else {
    sigma_eps_hat <- safe_mean("^sigma_eps\\[", expected_len = nt)
    L_corr_mean <- matrix(NA_real_, nt, nt)
    missing_L <- FALSE
    for (r in 1:nt) {
      for (c in 1:nt) {
        nm <- sprintf("L_corr[%d,%d]", r, c)
        if (nm %in% colnames(draws_df)) {
          L_corr_mean[r, c] <- mean(draws_df[[nm]])
        } else {
          missing_L <- TRUE
          L_corr_mean[r, c] <- NA_real_
        }
      }
    }
    if (missing_L || any(is.na(sigma_eps_hat))) {
      warning("Could not fully reconstruct Sigma: missing L_corr or sigma_eps elements.")
      Sigma_hat <- matrix(NA_real_, nt, nt)
    } else {
      L_Sigma_hat <- diag(sigma_eps_hat) %*% L_corr_mean
      Sigma_hat <- L_Sigma_hat %*% t(L_Sigma_hat)
    }
  }

  # Build population phi vector
  tau_own_hat   <- exp(tau_own_log_hat)
  tau_cross_hat <- exp(tau_cross_log_hat)
  lambda_own_hat_v   <- exp(lambda_own_log_hat)
  lambda_cross_hat_v <- exp(lambda_cross_log_hat)
  phi_pop_own <- tau_own_hat * lambda_own_hat_v * g_own_hat
  phi_pop_cross <- tau_cross_hat * lambda_cross_hat_v * g_cross_hat

  # Build own/cross indices
  idx_own <- integer(nt)
  idx_cross <- integer(nt * (nt - 1))
  {
    k_diag <- 1L
    k_cross <- 1L
    for (c in 1:nt) {
      for (r in 1:nt) {
        flat <- (c - 1) * nt + r
        if (r == c) {
          idx_own[k_diag] <- flat; k_diag <- k_diag + 1L
        } else {
          idx_cross[k_cross] <- flat; k_cross <- k_cross + 1L
        }
      }
    }
  }
  vec_phi_pop <- numeric(nt * nt)
  for (k in 1:nt) vec_phi_pop[idx_own[k]] <- phi_pop_own[k]
  for (k in 1:(nt * (nt - 1))) vec_phi_pop[idx_cross[k]] <- phi_pop_cross[k]

  # Subject-level phi0_j and Phi_j
  phi0_j_hat <- array(NA_real_, dim = c(J, nt))
  Phi_j_hat <- array(NA_real_, dim = c(J, nt, nt))
  for (j in 1:J) {
    z_phi0_j <- safe_mean(sprintf("^z_phi0\\[%d,", j), expected_len = nt)
    if (length(z_phi0_j) != nt) {
      warning("z_phi0[", j, "] has length ", length(z_phi0_j), " expected ", nt)
      z_phi0_j <- rep(NA_real_, nt)
    }
    phi0_j_hat[j, ] <- phi0_pop_hat + sigma_phi0_hat * z_phi0_j

    z_own_j <- safe_mean(sprintf("^z_own\\[%d,", j), expected_len = nt)
    if (length(z_own_j) != nt) {
      warning("z_own[", j, "] has length ", length(z_own_j), " expected ", nt)
      z_own_j <- rep(0, nt)
    }
    z_cross_j <- safe_mean(sprintf("^z_cross\\[%d,", j), expected_len = nt * (nt - 1))
    if (length(z_cross_j) != nt * (nt - 1)) {
      warning("z_cross[", j, "] has length ", length(z_cross_j), " expected ", nt * (nt - 1))
      z_cross_j <- rep(0, nt * (nt - 1))
    }

    delta_j <- numeric(nt * nt)
    for (k in 1:nt) delta_j[idx_own[k]] <- exp(sigma_re_own_log_hat) * z_own_j[k]
    for (k in 1:(nt * (nt - 1))) delta_j[idx_cross[k]] <- exp(sigma_re_cross_log_hat) * z_cross_j[k]

    vec_phi_j <- vec_phi_pop + delta_j
    mat_phi_j <- matrix(0, nt, nt)
    for (c in 1:nt) {
      for (r in 1:nt) {
        flat <- (c - 1) * nt + r
        mat_phi_j[r, c] <- vec_phi_j[flat]
      }
    }
    Phi_j_hat[j, , ] <- mat_phi_j
  }

  # Residuals from generated quantities: expect structure resid[j,d,t]
  resid_mean <- vector("list", J)
  any_resid_names <- any(grepl("^resid\\[", colnames(draws_df)))
  if (!any_resid_names) {
    warning("No resid[...] found in draws")
    for (j in 1:J) resid_mean[[j]] <- matrix(NA_real_, nrow = T, ncol = nt)
  } else {
    # Build J x nt x T mean array first
    resid_array <- array(NA_real_, dim = c(J, nt, T))
    for (j in 1:J) {
      for (d in 1:nt) {
        for (t in 1:T) {
          nm <- sprintf("resid[%d,%d,%d]", j, d, t)
          if (nm %in% colnames(draws_df)) {
            resid_array[j, d, t] <- mean(draws_df[[nm]])
          }
        }
      }
    }
    # Convert to list of T x nt per subject
    for (j in 1:J) {
      resid_mean[[j]] <- t(resid_array[j, , ])  # now T x nt
    }
  }

  list(
    phi0_pop = phi0_pop_hat,   # length nt
    phi0_j = phi0_j_hat,       # J x nt
    Phi_j = Phi_j_hat,         # J x nt x nt
    Sigma = Sigma_hat,         # nt x nt
    resid = resid_mean         # list of J matrices (T x nt)
  )
}
