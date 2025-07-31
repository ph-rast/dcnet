// file: stage2_garchhier.stan

data {
  int<lower=1> J;               // number of subjects
  int<lower=1> D;               // number of series per subject (D=nt)
  vector[J] chat[D];            // each subject’s fitted log‑c_h (intercept)
  vector[J] ahat[D];            // fitted a_h on original scale in (0,1)
  vector[J] bhat[D];            // fitted b_h on original scale in (0,1)
}

parameters {
  // log‑intercept hierarchies
  vector[D] mu_c;               // population mean of log‑c_h[d]
  vector<lower=0>[D] sigma_c;   // between‑subject SD of log‑c_h[d]

  // a_h on logit scale
  vector[D] mu_a_logit;         // pop mean on logit(a_h)
  vector<lower=0>[D] sigma_a;   // between‑subject SD on logit(a_h)

  // b_h on logit scale
  vector[D] mu_b_logit;         // pop mean on logit(b_h)
  vector<lower=0>[D] sigma_b;   // between‑subject SD on logit(b_h)
}

transformed parameters {
  // pull back onto (0,1) for a_h, b_h
  vector[D] mu_a = inv_logit(mu_a_logit);
  vector[D] mu_b = inv_logit(mu_b_logit);
}

model {
  // --- priors ---
  // c_h fixed (on log scale)
  mu_c    ~ normal(-0.9, 0.5);
  sigma_c ~ normal(0, 0.3);

  // a_h fixed  (on logit scale)   target mean ≈ 0.25
  mu_a_logit ~ normal(logit(0.25), 0.5);
  sigma_a    ~ normal(0, 0.3);

  // b_h fixed  (on logit scale)   target mean ≈ 0.37
  mu_b_logit ~ normal(logit(0.37), 0.5);
  sigma_b    ~ normal(0, 0.3);

  // --- likelihood ---
  for (d in 1:D) {
    chat[d] ~ normal(mu_c[d],    sigma_c[d]);
    ahat[d] ~ normal(mu_a[d],     sigma_a[d]);
    bhat[d] ~ normal(mu_b[d],     sigma_b[d]);
  }
}

generated quantities {
  // (optional) back‑transform the population means
  vector[D] mu_c_exp   = exp(mu_c);        // mean c_h on original scale
  vector[D] mu_a_raw   = mu_a;             // already in (0,1)
  vector[D] mu_b_raw   = mu_b;             // already in (0,1)
}
