
data {
  int<lower=1> J;                // subjects
  int<lower=1> T;                // time points
  int<lower=1> nt;               // number of series per subject
  array[J] matrix[T, nt] u;      // residuals from stage 1: T x nt per subject
}

parameters {
  // population-level fixed effects (on transformed scales)
  vector[nt] mu_c;               // log-variance intercepts
  vector[nt] mu_a_raw;           // raw for a before inv_logit
  vector[nt] mu_b_raw;           // raw for b before inv_logit (will be scaled to keep a+b < 0.99)

  // random effect standard deviations (hierarchy across subjects) 
  vector<lower=0>[nt] sd_c;
  vector<lower=0>[nt] sd_a;
  vector<lower=0>[nt] sd_b;

  // non-centered random effects
  array[J] vector[nt] z_c;       // for c_{j,d}
  array[J] vector[nt] z_a;       // for a_{j,d} raw
  array[J] vector[nt] z_b;       // for b_{j,d} raw
}

transformed parameters {
  matrix[J, nt] c_h;  // log-scale intercepts per subject-series
  matrix[J, nt] a_h;  // GARCH a coefficients
  matrix[J, nt] b_h;  // GARCH b coefficients

  for (j in 1:J) {
    for (d in 1:nt) {
      real c_raw = mu_c[d] + sd_c[d] * z_c[j][d];
      real a_raw = mu_a_raw[d] + sd_a[d] * z_a[j][d];
      real b_raw = mu_b_raw[d] + sd_b[d] * z_b[j][d];

      real a_val = inv_logit(a_raw);
      real b_val = (1 - a_val) * inv_logit(b_raw); 

      c_h[j,d] = c_raw;
      a_h[j,d] = a_val;
      b_h[j,d] = b_val;
    }
  }
}

model {
  // Priors
  mu_c ~ normal( -0.9, 1.0 );                 // prior on log-variance intercept (adjust as needed)
  mu_a_raw ~ normal( logit(0.2), 1.0 );        // prior centered near a ≈ 0.2
  mu_b_raw ~ normal( logit(0.3), 1.0 );        // prior centered near b ≈ 0.3

  sd_c ~ student_t(3, 0, 0.5);                 // weakly informative for variability across subjects
  sd_a ~ student_t(3, 0, 0.5);
  sd_b ~ student_t(3, 0, 0.5);

  // non-centered hierarchy
  for (j in 1:J) {
    to_vector(z_c[j]) ~ std_normal();
    to_vector(z_a[j]) ~ std_normal();
    to_vector(z_b[j]) ~ std_normal();
  }

  // Likelihood: GARCH recursion per subject j and series d
  // h2 is variance; initialize at unconditional variance for stationarity
  for (j in 1:J) {
    for (d in 1:nt) {
      vector[T] h2; // variance over time for series d of subject j

      // unconditional variance: exp(c) / (1 - a - b), safe because a+b<0.99
      real a_val = a_h[j,d];
      real b_val = b_h[j,d];
      real omega = exp(c_h[j,d]);
      real h2_uncond = omega / (1 - a_val - b_val);
      h2[1] = h2_uncond;

      for (t in 2:T) {
        real eps_prev_sq = square(u[j][t-1, d]);
        h2[t] = omega + a_val * eps_prev_sq + b_val * h2[t-1];
        // observation
        target += normal_lpdf(u[j][t, d] | 0, sqrt(h2[t]));
      }
      // optionally include first timepoint likelihood if you want to anchor:
      // target += normal_lpdf(u[j][1, d] | 0, sqrt(h2[1]));
    }
  }
}

generated quantities {
  matrix[J, nt] a_h_out = a_h;
  matrix[J, nt] b_h_out = b_h;
  matrix[J, nt] c_h_out = c_h;
}
