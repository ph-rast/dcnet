// STAGEn2: estimate between-subject SD of Fisher-z(S_j)

functions {
#include /functions/invvec.stan    // for generated quantities only
}

data {
  int<lower=1> J;             // number of subjects
  int<lower=2> nt;            // number of time series per subject
  int<lower=1> Sdim;          // = nt*(nt-1)/2  
  array[J] vector[Sdim] zhat;       // per-subject Fisher-z off-diag correlations
}

parameters {
  vector[Sdim] mu_z;          // population mean in z-space
  real<lower=0> sigma_z;      // between-subject SD
}

model {
  mu_z    ~ normal(0, 0.5);       // weakly informative around zero
  sigma_z ~ normal(0, 0.313);     // half-Normal mean approx. 0.25 if ranS_sd=0.25

  for (j in 1:J)
    zhat[j] ~ normal(mu_z, sigma_z);
}

generated quantities {
  corr_matrix[nt] R_pop;
  R_pop = invvec_to_corr( tanh(mu_z), nt );
}
