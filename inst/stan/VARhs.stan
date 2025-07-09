// DCC-Parameterization with random D components, random S fixed q's


data {
#include /data/data.stan
}

transformed data { 
  array[J] vector[nt] rts_m;
  array[J] vector[nt] rts_sd;
  int<lower=1> Sdim = (nt + nt*nt) %/% 2 - nt ; // Dimension of vec(S).
  //  int grainsize = 1;
  // Find positions of own and cross lags for flattened vec_phi
  array[nt]          int idx_own;
  array[nt * (nt-1)] int idx_cross;

  { // local scope to keep the counters private
    int k_diag  = 1;
    int k_cross = 1;

    // column-major sweep over the nt × nt matrix
    for (c in 1:nt) {
      for (r in 1:nt) {
        int flat = (c - 1) * nt + r;            // 1-based vector index
        if (r == c) {
          idx_own[k_diag]  = flat;
          k_diag          += 1;
        } else {
          idx_cross[k_cross] = flat;
          k_cross           += 1;
        }
      }
    }
  }

  
#include /transformed_data/xh_marker.stan
 
}
parameters {
  // VAR parameters
  
  // Horseshoe
  // semi-global scales
  real<lower=0> tau_own;
  real<lower=0> tau_cross;

  // local scales (Horseshoe example)
  vector<lower=0>[nt]   lambda_own;
  vector<lower=0>[nt*(nt-1)] lambda_cross;

  // standard-normal draws to be shrunk
  vector[nt]            g_own;
  vector[nt*(nt-1)]     g_cross;
  
  //#include /parameters/arma.stan
  // define arma.stan params right here:
  // BEGIN
  // Fixed effect vector
  vector[nt] phi0_fixed;

  cholesky_factor_corr[nt] phi0_L;
  vector<lower=0>[nt] phi0_tau; // ranef SD for phi0
  array[J] vector[nt] phi0_stdnorm; // Used to multiply with phi0_sd to obtain ranef on phi0
  // vec_phi_fixed is now defined in transformed params
  cholesky_factor_corr[nt*nt] phi_L;

  // Block separable SD's
  real<lower=0> sigma_re_own;
  real<lower=0> sigma_re_cross;

  // phi_tau now in transformed params
  // vector<lower=0>[nt*nt] phi_tau; // ranef SD for phi0

  array[J] vector[nt*nt] phi_stdnorm; // Used to multiply with phi0_sd to obtain ranef on phi0
  // END

  vector[nt] phi0_fixed2;  // In case of S_pred != NULL
  cholesky_factor_cov[nt] L_rescov;
  real< lower = 2 > nu; // nu for student_t 
}

transformed parameters {
  // transform vec_phi to nt*nt parameter matrix
  //matrix<lower = -1, upper = 1>[nt,nt] phi;
  array[J] vector[nt] phi0; // vector with fixed + random for intercept
  array[J] vector[nt*nt] phi; // vectorized VAR parameter matrix, fixed + random
  //  build vec_phi_fixed (length nt^2)
  // fixed effec: Grouped HS
  vector[nt*nt] vec_phi_fixed = rep_vector(0, nt*nt);
  cov_matrix[nt] rescov = multiply_lower_tri_self_transpose( L_rescov ) ;

  // Ranefs SD  vector
  vector<lower=0>[nt*nt] phi_tau = rep_vector(0, nt * nt);

  // Write own- and cross-lags into the right position in the flattened vectors:
  // own-lags (diagonal)
  for (k in 1:nt) {
    int pos = idx_own[k];                                  // 1-based flat index
    vec_phi_fixed[pos] = tau_own * lambda_own[k] * g_own[k];
    phi_tau[pos]       = sigma_re_own;
  }

  // cross-lags (off-diagonals)
  for (k in 1:(nt * (nt - 1))) {
    int pos = idx_cross[k];
    vec_phi_fixed[pos] = tau_cross * lambda_cross[k] * g_cross[k];
    phi_tau[pos]       = sigma_re_cross;
  }

  // --------- subject-specific coefficients -------------------------
  array[J,T] vector[nt] mu; 

  // VAR phi parameter

  // Initialize t=1
  for( j in 1:J){
    mu[j,1] = rts[j, 1]';
    for (t in 2:T){
      // Meanstructure model: DROP, if only VAR is allowed
#include /model_components/mu.stan 
    }
  }
}

model {
  // VAR
  phi0_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects
  phi_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects

  phi0_tau ~ normal(0, .1); // SD for multiplication with cholesky phi0_L
  phi_tau ~ normal(0, .1); // SD for multiplication with cholesky phi0_L

  // Horseshoe:
  //// standard normal prior
  g_own ~ std_normal();
  g_cross ~ std_normal();
  //// semi-global scales
  tau_own ~ cauchy(0, .5);
  tau_cross ~ cauchy(0, .25);
  /// local scales
  lambda_own   ~ normal(0, 1);
  lambda_cross ~ normal(0, 1);

  // VAR phi ranefs for own and cross lags:
  sigma_re_own   ~ normal(0, 0.1);   // looser
  sigma_re_cross ~ normal(0, 0.025);  // tight shrinkage

  // Prior for initial state
  // Cholesky factor prior (instead of Wishart)
  L_rescov ~ lkj_corr_cholesky(2);             // weak shrink toward I


  //if ( distribution == 1 )
  nu ~ normal( nt, 5 );
  //to_vector(phi0_fixed) ~ student_t(3,0,1);//multi_normal(rts_m, diag_matrix( rep_vector(1.0, nt) ) );
  phi0_fixed ~ student_t(3, 0, 1);
  vec_phi_fixed ~ student_t(3, 0, 1);
  
  // likelihood
  for( j in 1:J) {
    phi0_stdnorm[j] ~ std_normal();
    phi_stdnorm[j] ~ std_normal();
    if ( distribution == 0 ) {
      for(t in 1:T){
	to_vector(rts[j,t]) ~ multi_normal(mu[j,t], rescov);
      }
    } else if ( distribution == 1 ) {
      for(t in 1:T){ 
	to_vector(rts[j,t]) ~ multi_student_t(nu, mu[j,t], rescov); 
      } 
    } 
  }
}

generated quantities {
  //#include /generated/retrodict.stan
  array[J] matrix[T, nt] rts_out;
  array[J] vector[T] log_lik;


  if ( distribution == 0 ){
    for( j in 1:J ) {
      for (t in 1:T) {
	rts_out[j,t] = multi_normal_rng(mu[j,t], rescov)';
	log_lik[j,t] = multi_normal_lpdf(rts[j,t] | mu[j,t], rescov);
      }
    }
  } else if ( distribution == 1 ) {
    for( j in 1:J ) {
      for (t in 1:T) {
	rts_out[j,t] = multi_student_t_rng(nu, mu[j,t], rescov)';
	log_lik[j,t] = multi_student_t_lpdf(rts[j,t] | nu, mu[j,t], rescov);
      }
    }
  }
}

