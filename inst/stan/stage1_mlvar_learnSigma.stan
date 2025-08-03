// stage1_mlvar_learnSigma_rts.stan
data {
  int<lower=1> J;                              // subjects
  int<lower=1> T;                              // time length
  int<lower=1> nt;                             // series dimension
  real lbound;                                 // lower bound of observed
  real ubound;                                 // upper bound of observed
  array[J] matrix<lower = lbound, upper = ubound>[T, nt] rts; // T × nt per subject
}

transformed data {
  // own vs cross indices for vectorized nt×nt phi (column-major)
  array[nt] int idx_own;
  array[nt * (nt - 1)] int idx_cross;
  {
    int k_diag = 1;
    int k_cross = 1;
    for (c in 1:nt) {
      for (r in 1:nt) {
        int flat = (c - 1) * nt + r;
        if (r == c) {
          idx_own[k_diag] = flat;
          k_diag += 1;
        } else {
          idx_cross[k_cross] = flat;
          k_cross += 1;
        }
      }
    }
  }
}

parameters {
  // intercept hierarchy
  vector[nt] phi0_pop;
  vector[nt] sigma_phi0_log;
  array[J] vector[nt] z_phi0; // non-centered intercept deviations

  // population VAR phi with grouped horseshoe (own vs cross)
  real tau_own_log;
  real tau_cross_log;
  vector[nt]   lambda_own_log;
  vector[nt * (nt - 1)] lambda_cross_log;
  vector[nt]            g_own;
  vector[nt * (nt - 1)] g_cross;

  // subject-level deviations for phi: separate scales for own/cross
  real sigma_re_own_log;
  real sigma_re_cross_log;
  array[J] vector[nt]            z_own;
  array[J] vector[nt * (nt - 1)] z_cross;

  // residual covariance decomposition Sigma = D R D
  cholesky_factor_corr[nt] L_corr;    // correlation part
  vector<lower=0>[nt] sigma_eps;      // marginal scales
}

transformed parameters {
  vector[nt] sigma_phi0 = exp(sigma_phi0_log);
  // population phi
  real<lower=0> tau_own   = exp(tau_own_log);
  real<lower=0> tau_cross = exp(tau_cross_log);
  vector<lower=0>[nt]   lambda_own   = exp(lambda_own_log);
  vector<lower=0>[nt * (nt - 1)] lambda_cross = exp(lambda_cross_log);
  vector[nt] phi_pop_own   = tau_own .* lambda_own .* g_own;                    // own lags
  vector[nt * (nt - 1)] phi_pop_cross = tau_cross .* lambda_cross .* g_cross;   // cross lags

  vector[nt * nt] vec_phi_pop = rep_vector(0.0, nt * nt);
  for (k in 1:nt) {
    int pos = idx_own[k];
    vec_phi_pop[pos] = phi_pop_own[k];
  }
  for (k in 1:(nt * (nt - 1))) {
    int pos = idx_cross[k];
    vec_phi_pop[pos] = phi_pop_cross[k];
  }

  // subject-specific phi matrices and intercepts
  array[J] matrix[nt, nt] Phi;
  array[J] vector[nt] phi0_j;

  real<lower=0> sigma_re_own   = exp(sigma_re_own_log);
  real<lower=0> sigma_re_cross = exp(sigma_re_cross_log);

  for (j in 1:J) {
    // intercept non-centered
    phi0_j[j] = phi0_pop + sigma_phi0 .* z_phi0[j];

    // build subject delta
    vector[nt * nt] delta_j = rep_vector(0.0, nt * nt);
    for (k in 1:nt) {
      int pos = idx_own[k];
      delta_j[pos] = sigma_re_own * z_own[j][k];
    }
    for (k in 1:(nt * (nt - 1))) {
      int pos = idx_cross[k];
      delta_j[pos] = sigma_re_cross * z_cross[j][k];
    }

    vector[nt * nt] vec_phi_j = vec_phi_pop + delta_j;

    // reshape into matrix (column-major)
    for (c in 1:nt) {
      for (r in 1:nt) {
        int flat = (c - 1) * nt + r;
        Phi[j][r, c] = vec_phi_j[flat];
      }
    }
  }

  // build Sigma
  matrix[nt, nt] L_Sigma = diag_pre_multiply(sigma_eps, L_corr);
  matrix[nt, nt] Sigma = multiply_lower_tri_self_transpose(L_Sigma); // residual covariance
}

model {
  // ---- priors ----
  // intercept hierarchy
  phi0_pop ~ normal(0, 1);
  sigma_phi0 ~ normal(0, 0.5);
  for (j in 1:J)
    to_vector(z_phi0[j]) ~ std_normal();

  // horseshoe-style shrinkage for population phi
  tau_own_log   ~ student_t(3, log(0.3), 0.5);
  tau_cross_log ~ student_t(3, log(0.05), 0.2);
  lambda_own_log   ~ student_t(3, 0, 1);
  lambda_cross_log ~ student_t(3, 0, 1);
  to_vector(g_own)   ~ std_normal();
  to_vector(g_cross) ~ std_normal();

  // subject deviations
  sigma_re_own_log   ~ normal(log(0.01), 0.5);
  sigma_re_cross_log ~ normal(log(0.01), 0.5);
  for (j in 1:J) {
    to_vector(z_own[j]) ~ std_normal();
    to_vector(z_cross[j]) ~ std_normal();
  }

  // Sigma priors: regularize residual covariance
  L_corr ~ lkj_corr_cholesky(2);
  sigma_eps ~ student_t(3, 0, 1);

  // ---- likelihood: VAR(1) ----
  for (j in 1:J) {
    for (t in 2:T) {
      vector[nt] rts_lag = to_vector(rts[j][t - 1]); // previous time point as column
      vector[nt] mu = Phi[j] * (rts_lag - phi0_j[j])  + phi0_j[j]; // centered version with rts_lag - phi0
      vector[nt] rts_t = to_vector(rts[j][t]);
      rts_t ~ multi_normal(mu, Sigma);
    }
  }
}

generated quantities {
  array[J] matrix[nt, T] resid;  // residuals: rts_{j,t} - [phi0_j + Phi_j * (rts_{j,t-1} - phi0_j)]
 
  // Residuals (centered)
  for (j in 1:J) {
    vector[nt] phi0_j_loc = phi0_pop + sigma_phi0 .* z_phi0[j]; // reuse same transform
    for (t in 2:T) {
      vector[nt] rts_lag = to_vector(row(rts[j], t - 1));
      vector[nt] mu = phi0_j_loc + Phi[j] * (rts_lag - phi0_j_loc);
      vector[nt] rts_t = to_vector(row(rts[j], t));
      resid[j][, t] = rts_t - mu;
    }
    // first timepoint residuals zero
    for (d in 1:nt) {
      resid[j][d, 1] = 0;
    }
  }
}
