// DCC-Parameterization with random D components, random S fixed q's
functions {
#include /functions/jacobian.stan
#include /functions/invvec.stan
#include /functions/cov2cor.stan 

}

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
}


parameters {
  // VAR parameters
  // Horseshoe
  // semi-global scales
  real tau_own_log;
  real tau_cross_log;

  // local scales (Horseshoe)
  vector[nt]   lambda_own_log;
  vector[nt*(nt-1)] lambda_cross_log;

  // standard-normal draws to be shrunk
  vector[nt]            g_own;
  vector[nt*(nt-1)]     g_cross;
  
  //#include /parameters/arma.stan
  // define arma.stan params right here:
  // BEGIN
  // Fixed effect vector
  vector[nt] phi0_fixed;

  cholesky_factor_corr[nt] phi0_L;
  vector[nt] phi0_tau_log; // ranef SD for phi0
  array[J] vector[nt] phi0_stdnorm; // Used to multiply with phi0_sd to obtain ranef on phi0
  // vec_phi_fixed is now defined in transformed params
  cholesky_factor_corr[nt*nt] phi_L;

  // Block separable SD's
  real sigma_re_own_log;
  real sigma_re_cross_log;

  // phi_tau now in transformed params
  // vector<lower=0>[nt*nt] phi_tau; // ranef SD for phi0

  array[J] vector[nt*nt] phi_stdnorm; // Used to multiply with phi0_sd to obtain ranef on phi0
  // END

  vector[nt] phi0_fixed2;  // In case of S_pred != NULL 
  
  // predictor for H
#include /parameters/predH.stan

  // GARCH h parameters on variance metric
  vector[nt] c_h_fixed; // variance on log metric
  vector[nt] a_h_fixed;
  vector[nt] b_h_fixed;

  // Random effects for c_h
  cholesky_factor_corr[nt] c_h_L;
  vector[nt] c_h_tau_log;
  array[J] vector[nt] c_h_stdnorm;
  
  // Random effects fo a_h
  cholesky_factor_corr[nt] a_h_L;
  vector[nt] a_h_tau_log; // ranef SD for phi0
  array[J] vector[nt] a_h_stdnorm; // Used to multiply with phi0_sd to obtain ranef on phi0

  // Random effects fo b_h
  cholesky_factor_corr[nt] b_h_L;
  vector[nt] b_h_tau_log; // ranef SD for phi0
  array[J] vector[nt] b_h_stdnorm; // Used to multiply with phi0_sd to obtain ranef on phi0

  // Create raws to facilitate prior assignments
  vector<lower=0, upper=1>[nt] a_h_pop_raw;
  vector<lower=0, upper=1>[nt] b_h_pop_raw_raw;

  
  // GARCH q parameters
  //real l_a_q; // Define on log scale so that it can go below 0
  
  real l_a_q_sigma_log; //random effect variance
  array[J] real l_a_q_stdnorm;
  
  //real<lower=0, upper = (1 - a_q) > b_q; //
  //real l_b_q; 
  real l_b_q_sigma_log; //random effect variance
  array[J] real l_b_q_stdnorm;

  // Try with a_q/b_q params that are on original scale for priors:
  real<lower=0, upper=1> a_q_pop_raw;         // population mean of a_q
  real<lower=0, upper=1> b_q_pop_raw_raw;     // temporary, to enforce sum<1

  
  corr_matrix[nt] S;
  
  corr_matrix[nt] S2;

  // R1 init
  array[J] corr_matrix[nt] R1_init;
  
  // Qr1 init
  array[J] cov_matrix[nt] Qr1_init;
  // D1 init
  array[J] vector<lower = 0>[nt] D1_init;
  // u1 init
  array[J] vector[nt] u1_init;

  real< lower = 2 > nu; // nu for student_t
}

transformed parameters {
  real<lower=0> sigma_re_own = exp(sigma_re_own_log);
  real<lower=0> sigma_re_cross = exp(sigma_re_cross_log);
  real<lower=0> tau_own = exp(tau_own_log);
  real<lower=0> tau_cross = exp(tau_cross_log);
  vector<lower=0>[nt]   lambda_own = exp(lambda_own_log);
  vector<lower=0>[nt*(nt-1)] lambda_cross = exp(lambda_cross_log);
  vector<lower=0>[nt] phi0_tau = exp(phi0_tau_log); // ranef SD for phi0
  vector<lower=0>[nt] c_h_tau = exp(c_h_tau_log);
  vector<lower=0>[nt] a_h_tau = exp(a_h_tau_log); // ranef SD for phi0
  vector<lower=0>[nt] b_h_tau = exp(b_h_tau_log); // ranef SD for phi0
  real<lower=0> l_a_q_sigma = exp(l_a_q_sigma_log); //random effect variance
  real<lower=0> l_b_q_sigma = exp(l_b_q_sigma_log); //random effect variance
  // transform vec_phi to nt*nt parameter matrix
  //matrix<lower = -1, upper = 1>[nt,nt] phi;
  array[J] vector[nt] phi0; // vector with fixed + random for intercept
  array[J] vector[nt*nt] phi; // vectorized VAR parameter matrix, fixed + random
  //  build vec_phi_fixed (length nt^2)
  // fixed effec: Grouped HS
  vector[nt*nt] vec_phi_fixed = rep_vector(0, nt*nt);

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

  
  // shrink b_q so that a_q + b_q ≤ 0.99   (add a small error margin)
  real<lower=0, upper=1> a_q_pop = a_q_pop_raw;
  real<lower=0, upper=1> b_q_pop = (1 - a_q_pop - 0.01) * b_q_pop_raw_raw;

  vector<lower=0, upper=1>[nt] a_h_pop;
  vector<lower=0, upper=1>[nt] b_h_pop;
  for (d in 1:nt) {
    a_h_pop[d] = 0.50 * a_h_pop_raw[d];                        // centred ~0.25
    b_h_pop[d] = (0.99 - a_h_pop[d]) * b_h_pop_raw_raw[d];     // centred ~0.37
  }
  // --------- subject-specific coefficients -------------------------
  array[J] vector<lower=0, upper=1>[nt] a_h_sum;
  array[J] vector<lower=0, upper=1>[nt] b_h_sum;
  
  for (j in 1:J) {
    vector[nt] a_re = a_h_tau .* a_h_stdnorm[j];               // non-centred RE
    vector[nt] b_re = b_h_tau .* b_h_stdnorm[j];
    
    // add REs on logit scale, then transform back  
    a_h_sum[j] =
      inv_logit( logit(a_h_pop) + a_re );                      // stays in (0,1)
    
    // ensure a_h + b_h < 0.99 element-wise
    for (d in 1:nt) {
      real b_raw = inv_logit( logit(b_h_pop[d]) + b_re[d] );   // (0,1)
      b_h_sum[j,d] = (0.99 - a_h_sum[j,d]) * b_raw;            // (0, 0.99-a)
    }
  }


  

  
  // logit versions used elsewhere in my code
  real l_a_q = logit(a_q_pop);
  real l_b_q = logit(b_q_pop);
  
  //Scale
  array[J] real<lower = 0, upper = 1> a_q;
  array[J] real<lower = 0, upper = 1> b_q;

  array[J] vector[nt] c_h;
  array[J] vector[nt] c_h_random; // variance on log metric
  //vector[nt] a_h[J];
  array[J] vector[nt] a_h_random; // variance on log metric
  array[J] vector[nt] b_h_random;

  //array[J] vector<lower=0, upper = 1>[nt] a_h_sum;
  //array[J] vector<lower=0, upper = 1>[nt] b_h_sum;
  
  array[J,T] cov_matrix[nt] H;
  array[J,T] corr_matrix[nt] R;
  array[J,T-1] vector[nt] rr;
  array[J,T] vector[nt] mu;
  array[J,T] vector[nt] D;
  array[J,T] cov_matrix[nt] Qr;
  array[J,T] vector[nt] Qr_sdi;
  array[J,T] vector[nt] u;
  
  array[J] vector<lower = 0>[nt] vd;
  array[J] vector<lower = 0>[nt] ma_d;
  array[J] vector<lower = 0>[nt] ar_d;
  

  array[J] real l_a_q_r;
  array[J] real l_b_q_r;
  
  // VAR phi parameter

  // Initialize t=1
  for( j in 1:J){
    
    mu[j,1] = rts[j, 1]'; //mu[,1] is not used, just copy rts1 as starting value 
    D[j,1] = D1_init[j];//sqrt( diagonal((rts[j]' * rts[j]) / (nt - 1)) );//D1_init[j];//
    u[j,1] = u1_init[j];
    Qr[j,1] = Qr1_init[j];//

    //Qr[j,1] = ((rts[j]' * rts[j]) / (nt - 1));// Werden alle flach
    Qr_sdi[j,1] = 1 ./ sqrt(diagonal(Qr[j,1])); //
    
    //R[j,1] = R1_init[j];//cov2cor((rts[j]' * rts[j]) / (nt - 1));
    //R[j,1]=diag_matrix(rep_vector(1.0, nt));
    R[j,1] = quad_form_diag(Qr[j,1], Qr_sdi[j,1]); //
    H[j,1] = quad_form_diag(R[j,1],  D[j,1]);

    ////////////
    // Ranefs //
    ////////////
    
    // R part
    l_a_q_r[j] = l_a_q_sigma * l_a_q_stdnorm[j];
    l_b_q_r[j] = l_b_q_sigma * l_b_q_stdnorm[j];
    //a_q[j] =          1 ./ ( 1 + exp(-(l_a_q + l_a_q_r[j])) );
    // Replaced by stans inv_logit function:
    a_q[j] =          inv_logit( l_a_q + l_a_q_r[j] );
    //b_q[j] = (1-a_q[j]) ./ ( 1 + exp(-(l_b_q + l_b_q_r[j])) );
    b_q[j] = (1 - a_q[j]) * inv_logit(l_b_q + l_b_q_r[j]);
    // D part
    // Should random effect corrleations be estimated:
    // Full estimation of corrmat with simplify_ch == 0
    // Estimate only diagonal with simplify_ch == 1
    if(simplify_ch==0){
      c_h_random[j] = (diag_pre_multiply(c_h_tau, c_h_L)*c_h_stdnorm[j]);
    } else if(simplify_ch == 1) {
      c_h_random[j] = c_h_tau .* c_h_stdnorm[j];
    }
    c_h[j] = c_h_fixed + c_h_random[j];

    // Same approach here with simplify_ah
    if(simplify_ah==0){
      a_h_random[j] = diag_pre_multiply(a_h_tau, a_h_L)*a_h_stdnorm[j];
    } else if(simplify_ah == 1) {
      a_h_random[j] = a_h_tau .* a_h_stdnorm[j];
    }
    
    // Bound sum of fixed and ranef between 0 and 1 with logistic function
    //a_h_sum[j] = 1 ./ (1 + exp(-( a_h_fixed + a_h_random[j] )) );
    a_h_sum[j] = inv_logit(a_h_fixed + a_h_random[j]);
    
    if(simplify_bh==0){      
      b_h_random[j] = (diag_pre_multiply(b_h_tau, b_h_L)*b_h_stdnorm[j]);
    } else if(simplify_bh == 1) {
      b_h_random[j] = b_h_tau .* b_h_stdnorm[j];
    }
    // Bound sum of fixed and ranef between 0 and 1
     b_h_sum[j] = (1 - a_h_sum[j]) ./ (1 + exp(-(b_h_fixed + b_h_random[j])) );
    //b_h_sum[j] = (1 - a_h_sum[j]) .* inv_logit(b_h_fixed + b_h_random[j]);
    /////////////////////
    // Model Structure //
    /////////////////////
    
    for (t in 2:T){
      // Meanstructure model: DROP, if only VAR is allowed
#include /model_components/mu.stan
      
      
      for(d in 1:nt){
	vd[j,d]   = 0.0;
	ma_d[j,d] = 0.0;
	ar_d[j,d] = 0.0;
	// GARCH MA component
	for (q in 1:min( t-1, Q) ) {
	  rr[j, t-q, d] = square( rts[j, t-q, d] - mu[j, t-q, d] );
	  ma_d[j, d] = ma_d[j, d] + a_h_sum[j, d]*rr[j, t-q, d] ;
	}
	// print("ma_d: ", "TS:", d, " Value:", ma_d[d], " T:", t);
	// GARCH AR component
	for (p in 1:min( t-1, P) ) {
	  ar_d[j,d] = ar_d[j,d] + b_h_sum[j, d]*D[j, t-p, d]^2;
	}
	// print("ar_d: ", "TS:", d, " Value:", ar_d[d], " T:", t);
	// Predictor on diag (given in xC)
	//if ( xC_marker >= 1) {
	//  vd[j,d] = exp( c_h[j,d] + beta[d] * xC[t, d] ) + ma_d[j,d] + ar_d[j,d];
	//} else if ( xC_marker == 0) {
	  vd[j,d] = exp( c_h[j,d] )  + ma_d[j,d] + ar_d[j,d];
	  //}

	D[j, t, d] = sqrt( vd[j,d] );
      }
      u[j,t,] = diag_matrix(D[j,t]) \ (rts[j,t]'- mu[j,t]) ; // cf. comment about taking inverses in stan manual p. 482 re:Inverses - inv(D)*y = D \ a

      //Introduce predictor for S (time-varying)
      if (S_pred[j,t] == 0){
	Qr[j,t ] = (1 - a_q[j] - b_q[j]) * S + a_q[j] * (u[j, t-1 ] * u[j, t-1 ]') + b_q[j] * Qr[j, t-1]; // S and UU' define dimension of Qr
      } else if (S_pred[j,t] == 1){
	Qr[j,t ] = (1 - a_q[j] - b_q[j]) * S2 + a_q[j] * (u[j, t-1 ] * u[j, t-1 ]') + b_q[j] * Qr[j, t-1];
      }
      Qr_sdi[j, t] = 1 ./ sqrt(diagonal(Qr[j, t])) ; // inverse of diagonal matrix of sd's of Qr
      //    R[t,] = quad_form_diag(Qr[t,], inv(sqrt(diagonal(Qr[t,]))) ); // Qr_sdi[t,] * Qr[t,] * Qr_sdi[t,];
      R[j,t] = quad_form_diag(Qr[j,t], Qr_sdi[j,t]); // 
      H[j,t] = quad_form_diag(R[j,t],  D[j,t]);  // H = DRD;
    }
  }
}

model {
  // print("Upper Limits:", UPs);
  // population-level priors
  a_q_pop_raw     ~ beta(1.8, 8.2);   // mean ≈ 0.09,  95 % ≈ (0.01, 0.23)
  b_q_pop_raw_raw ~ beta(3.8,  6.2);   // mean ≈ 0.25 BEFORE scaling
  // population means: weak Beta priors (most mass ≤ 0.3)
  a_h_pop_raw      ~ beta(1.5, 18.5);         // mean ≈ 0.20
  b_h_pop_raw_raw  ~ beta(3.4, 16.6);         // mean ≈ 0.25 (scaled later)
  
  
  // priors
  //l_a_q ~ student_t(3, -1.5, 2);
  //l_b_q ~ student_t(3, -1.5, 2);
  l_a_q_sigma_log ~ normal(log(0.63), 0.5); //student_t(3, 0, 2);
  to_vector(l_a_q_stdnorm) ~ std_normal();
  l_b_q_sigma_log ~ normal(log(0.63), .5); //student_t(3, 0, 2);
  to_vector(l_b_q_stdnorm) ~ std_normal();

  // VAR
  phi0_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects
  phi_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects
  //D part in DRD
  c_h_L ~ lkj_corr_cholesky(1);
  a_h_L ~ lkj_corr_cholesky(1);
  b_h_L ~ lkj_corr_cholesky(1);
  // R part in DRD

  phi0_tau_log ~ normal(log(0.5), .5); // SD for multiplication with cholesky phi0_L
  //phi_tau ~ normal(0, 0.0063); // SD for multiplication with cholesky phi0_L
  c_h_tau_log ~ normal(-1, .5); // SD for c_h ranefs
  a_h_tau_log ~ normal(-2, .5); // SD for c_h ranefs
  b_h_tau_log ~ normal(-2, .5);

  // Horseshoe:
  //// standard normal prior
  g_own ~ std_normal();
  g_cross ~ std_normal();
  //// semi-global scales
  tau_own_log ~ student_t(3, log(0.3), .5); //cauchy(0, 0.5);
  tau_cross_log ~ student_t(3, log(0.05), .2); // cauchy(0, 0.02);
  /// local scales
  lambda_own_log   ~ student_t(3, 0, 1);
  lambda_cross_log ~ student_t(3, 0, 1);

  // VAR phi ranefs for own and cross lags:
  sigma_re_own_log   ~ normal(log(0.01), 0.5);  
  sigma_re_cross_log ~ normal(log(0.01), 0.5);  
  
  // C
  to_vector(beta) ~ std_normal();
  to_vector(c_h_fixed) ~ normal( -0.9, .2);
  to_vector(a_h_fixed) ~ normal( -2.5, .2);
  to_vector(b_h_fixed) ~ normal( -1.5, .2);
  // Prior for initial state
  
  // Prior on nu for student_t
  //if ( distribution == 1 )
  nu ~ normal( nt, 5 );
  //to_vector(theta) ~ std_normal();
  //to_vector(phi) ~ std_normal();
  phi0_fixed ~ student_t(3, 0, 1);
  // vec_phi_fixed ~ normal(0, 5);
  //  to_vector(a_h) ~ normal(0, .5);
  //to_vector(b_h) ~ normal(0, .5);
  //S ~ lkj_corr( 1 );
  S2 ~ lkj_corr( 1 );
  S ~ lkj_corr( 1 );
  
  // likelihood
  for( j in 1:J) {
    R1_init[j] ~ lkj_corr( 1 );
    to_vector(D1_init[j]) ~ lognormal(-1, 1);
    to_vector(u1_init[j]) ~ std_normal();
    Qr1_init[j] ~ wishart(nt + 1.0, diag_matrix(rep_vector(1.0, nt)) );
    // UL transform jacobian
    /* for(k in 1:nt) { */
    /*   ULs[j,k] ~ uniform(0, UPs[j,k]); // Truncation not needed. */
    /*   target += a_b_scale_jacobian(0.0, ULs[j,k], b_h_sum_s[j,k]); */
    /* } */
    phi0_stdnorm[j] ~ std_normal();
    phi_stdnorm[j] ~ std_normal();
    c_h_stdnorm[j] ~ std_normal();
    a_h_stdnorm[j] ~ std_normal();
    b_h_stdnorm[j] ~ std_normal();
      // Likelihood
    if ( distribution == 0 ) {
      for(t in 1:T){
	rts[j,t] ~ multi_normal(mu[j,t], H[j,t]);
      }
    } else if ( distribution == 1 ) {
      for(t in 1:T){
	rts[j,t] ~ multi_student_t(nu, mu[j,t], H[j,t]);
      }
    }
  }
}

generated quantities {
  /* array[J,T] cov_matrix[nt] precision; */
  /* array[J,T] matrix[nt,nt] pcor; */
  /* matrix[nt,nt] Sfixed; */
#include /generated/retrodict.stan
  /* for(j in 1:J){ */
  /*   for(t in 1:T){ */
  /*     precision[j,t] = inverse(H[j,t]); */
  /*     pcor[j, t] = - cov2cor(precision[j,t]); */
  /*   } */
  /* } */
}
