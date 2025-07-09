// DCC-Parameterization with random D components, random S fixed q's
functions { 
  //#include /functions/jacobian.stan
#include /functions/invvec.stan
#include /functions/cov2cor.stan
// reduced sums:  
  real partial_log_lik(array[] matrix rts_slice,
                       int  start, int end,
                       array[,] vector mu,
                       array[,] matrix L_H,
                       int  distribution,
                       real nu) {

    int T = dims(rts_slice[1])[1];        // # time points
    real lp = 0;

    for (i in 1:size(rts_slice)) {
      int j = start + i - 1;              // global subject index

      for (t in 1:T) {
        vector[num_elements(rts_slice[i][1])] y =
            to_vector(rts_slice[i][t]');  // row to column vector

        if (distribution == 0) {
          lp += multi_normal_cholesky_lpdf(y | mu[j, t], L_H[j, t]);
        } else {
          lp += multi_student_t_cholesky_lpdf(y | nu, mu[j, t], L_H[j, t]);
        }
      }
    }
    return lp;
  }
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

  
#include /transformed_data/xh_marker.stan
 
  /* if( meanstructure == 0 ){ */
  /*   for ( j in 1:J ){ */
  /*     for ( i in 1:nt ){ */
  /* 	rts_m[j] = rep_vector(mean(rts[j,i]),nt); */
  /* 	rts_sd[j] =  rep_vector(sd(rts[j,i]),nt); */
  /*     } */
  /*   } */
  /* } else if (meanstructure == 1 || meanstructure == 2 ){ */
  /*   // set rts_m to first element in ts */
  /*   for ( j in 1:J ){ */
  /*     rts_m[j] = rts[j,1]'; */
  /*     rts_sd[j] = rep_vector(sd(rts[1,1]),nt); */
  /*   } */
  /* } */
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
  
  // predictor for H
#include /parameters/predH.stan

  // GARCH h parameters on variance metric
  vector[nt] c_h_fixed; // variance on log metric
  vector[nt] a_h_fixed;
  vector[nt] b_h_fixed;

  // Random effects for c_h
  cholesky_factor_corr[nt] c_h_L;
  vector<lower=0>[nt] c_h_tau;
  array[J] vector[nt] c_h_stdnorm;
  
  // Random effects fo a_h
  cholesky_factor_corr[nt] a_h_L;
  vector<lower=0>[nt] a_h_tau; // ranef SD for phi0
  array[J] vector[nt] a_h_stdnorm; // Used to multiply with phi0_sd to obtain ranef on phi0

  // Random effects fo b_h
  cholesky_factor_corr[nt] b_h_L;
  vector<lower=0>[nt] b_h_tau; // ranef SD for phi0
  array[J] vector[nt] b_h_stdnorm; // Used to multiply with phi0_sd to obtain ranef on phi0

  // Create raws to facilitate prior assignments
  vector<lower=0, upper=1>[nt] a_h_pop_raw;
  vector<lower=0, upper=1>[nt] b_h_pop_raw_raw;

  
  // GARCH q parameters
  //real l_a_q; // Define on log scale so that it can go below 0
  
  real<lower=0> l_a_q_sigma; //random effect variance
  array[J] real l_a_q_stdnorm;
  
  //real<lower=0, upper = (1 - a_q) > b_q; //
  //real l_b_q; 
  real<lower=0> l_b_q_sigma; //random effect variance
  array[J] real l_b_q_stdnorm;

  // Try with a_q/b_q params that are on original scale for priors:
  real<lower=0, upper=1> a_q_pop_raw;         // population mean of a_q
  real<lower=0, upper=1> b_q_pop_raw_raw;     // temporary, to enforce sum<1

  
  //corr_matrix[nt] S;
  array[J] vector[Sdim] S_vec_stdnorm; 
  vector<lower=0>[Sdim] S_vec_tau; 
  vector[Sdim] S_vec_fixed;  // Vectorized fixed effect for S
  vector[Sdim] S_vec_fixed2; // Vectorized fixed effect for S
  array[J] vector[Sdim] S_vec_sd;
  
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

  
  // shrink b_q so that a_q + b_q ≤ 0.99   (add a small ε margin)
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

  //old
  //array[J,T] cov_matrix[nt] H;
  // try cholesky
  array[J,T] cholesky_factor_cov[nt] L_H;
  
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
  

  // fixed + ranef vector for vec(S) part
  array[J] vector[Sdim] S_Lv;
  array[J] corr_matrix[nt] S;

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
    //H[j,1] = quad_form_diag(R[j,1],  D[j,1]);
    matrix[nt,nt] Hi1 = quad_form_diag(R[j,1], D[j,1]);
    Hi1 = 0.5 * (Hi1 + Hi1');
    Hi1 += diag_matrix(rep_vector(1e-6, nt));
    L_H[j,1] = cholesky_decompose(Hi1);
    
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
	if ( xC_marker >= 1) {
	  vd[j,d] = exp( c_h[j,d] + beta[d] * xC[t, d] ) + ma_d[j,d] + ar_d[j,d];
	} else if ( xC_marker == 0) {
	  vd[j,d] = exp( c_h[j,d] )  + ma_d[j,d] + ar_d[j,d];
	}

	D[j, t, d] = sqrt( vd[j,d] );
      }
      u[j,t,] = diag_matrix(D[j,t]) \ (rts[j,t]'- mu[j,t]) ; // cf. comment about taking inverses in stan manual p. 482 re:Inverses - inv(D)*y = D \ a
      //u[j, t, ] = (rts[j, t]' - mu[j, t]) ./ D[j, t];

      // Second S_lv stuff is just a placeholder: TODO, add S_Lv2[j] or drop completely
      if (S_pred[j,t] == 0){
	S_Lv[j] = tanh( S_vec_fixed  + S_vec_tau .*S_vec_stdnorm[j] );
      } else if (S_pred[j,t] == 1){
	S_Lv[j] = tanh( S_vec_fixed2 + S_vec_tau .*S_vec_stdnorm[j] );
      }
      S[j] = invvec_to_corr(S_Lv[j], nt);

      //Introduce predictor for S (time-varying)
      Qr[j,t ] = (1 - a_q[j] - b_q[j]) * S[j] + a_q[j] * (u[j, t-1 ] * u[j, t-1 ]') + b_q[j] * Qr[j, t-1]; // S and UU' define dimension of Qr
      Qr_sdi[j, t] = 1 ./ sqrt(diagonal(Qr[j, t])) ; // inverse of diagonal matrix of sd's of Qr
      R[j,t] = quad_form_diag(Qr[j,t], Qr_sdi[j,t]); // 
      // H[j,t] = quad_form_diag(R[j,t],  D[j,t]);  // H = DRD;

      /* Add Jitter to avoid non SPD matrix error: */
      matrix[nt,nt] Hij = quad_form_diag(R[j,t],  D[j,t]);
      // enforce symmetry
      Hij = 0.5 * (Hij + Hij');
      
      // add 1e-6 to diagonal to guarantee SPD
      Hij +=  diag_matrix(rep_vector(1e-6, nt));
      L_H[j,t] = cholesky_decompose(Hij);
      
      
      //H[j,t] = 0.5 * (H[j,t] + H[j,t]');
    }
  }
}

model {
  // print("Upper Limits:", UPs);
  // population-level priors
  a_q_pop_raw     ~ beta(2, 20);   // mean ≈ 0.09,  95 % ≈ (0.01, 0.23)
  b_q_pop_raw_raw ~ beta(2,  6);   // mean ≈ 0.25 BEFORE scaling
  // population means: weak Beta priors (most mass ≤ 0.3)
  a_h_pop_raw      ~ beta(2, 8);         // mean ≈ 0.20
  b_h_pop_raw_raw  ~ beta(2, 6);         // mean ≈ 0.25 (scaled later)
  
  
  // priors
  //l_a_q ~ student_t(3, -1.5, 2);
  //l_b_q ~ student_t(3, -1.5, 2);
  l_a_q_sigma ~ normal(0, 0.05); //student_t(3, 0, 2);
  to_vector(l_a_q_stdnorm) ~ std_normal();
  l_b_q_sigma ~ normal(0, 0.05); //student_t(3, 0, 2);
  to_vector(l_b_q_stdnorm) ~ std_normal();

  // VAR
  phi0_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects
  phi_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects
  //D part in DRD
  c_h_L ~ lkj_corr_cholesky(1);
  a_h_L ~ lkj_corr_cholesky(1);
  b_h_L ~ lkj_corr_cholesky(1);
  // R part in DRD

  phi0_tau ~ normal(0, .1); // SD for multiplication with cholesky phi0_L
  phi_tau ~ normal(0, .1); // SD for multiplication with cholesky phi0_L
  c_h_tau ~ normal(0, .05); // SD for c_h ranefs
  a_h_tau ~ normal(0, .05); // SD for c_h ranefs
  b_h_tau ~ normal(0, .05);

  // Horseshoe:
  //// standard normal prior
  g_own ~ std_normal();
  g_cross ~ std_normal();
  //// semi-global scales
  tau_own ~ cauchy(0, 0.5);
  tau_cross ~ cauchy(0, 0.02);
  /// local scales
  lambda_own   ~ normal(0, 1);
  lambda_cross ~ normal(0, 1);

  // VAR phi ranefs for own and cross lags:
  sigma_re_own   ~ normal(0, 0.1);   // looser
  sigma_re_cross ~ normal(0, 0.025);  // tight shrinkage
  
  // C
  to_vector(beta) ~ std_normal();
  to_vector(c_h_fixed) ~ normal( -2, .11);
  to_vector(a_h_fixed) ~ normal( -2, .11);
  to_vector(b_h_fixed) ~ normal( -2, .11);
  // Prior for initial state
  
  // Prior on nu for student_t
  //if ( distribution == 1 )
  nu ~ normal( nt, 5 );
  //to_vector(phi0_fixed) ~ student_t(3,0,1);//multi_normal(rts_m, diag_matrix( rep_vector(1.0, nt) ) );
  phi0_fixed ~ student_t(3, 0, 1);
  vec_phi_fixed ~ student_t(3, 0, 1);
  S_vec_fixed ~ std_normal();
  S_vec_fixed2 ~ std_normal();  
  
  // likelihood
  for( j in 1:J) {
    S_vec_tau ~ student_t(3, 0, 3);
    S_vec_stdnorm[j] ~ std_normal();
  
    //R1_init[j] ~ lkj_corr( 1 );
    to_vector(D1_init[j]) ~ lognormal(c_h[j], 1);
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
      target += reduce_sum(partial_log_lik,
			   rts,
			   grainsize,
			   mu, L_H,
			   distribution, nu);
  }
}

generated quantities {
  //array[J,T] cov_matrix[nt] precision;
  array[J,T] matrix[nt,nt] pcor;
  matrix[nt,nt] Sfixed;
  matrix[nt,nt] Sfixed2;
  array[J,T] cov_matrix[nt] H; 
  //#include /generated/retrodict.stan
  array[J] matrix[T, nt] rts_out;
  array[J] vector[T] log_lik;

  if ( distribution == 0 ){
    for( j in 1:J ) {
      for (t in 1:T) {
	H[j, t] = multiply_lower_tri_self_transpose(L_H[j,t]);
	rts_out[j,t] = multi_normal_rng(mu[j,t],
					multiply_lower_tri_self_transpose(L_H[j,t]))';
	log_lik[j,t] = multi_normal_lpdf(rts[j,t] | mu[j,t],
					 multiply_lower_tri_self_transpose(L_H[j,t]));
      }
    }
  } else if ( distribution == 1 ) {
    for( j in 1:J ) {
      for (t in 1:T) {
	rts_out[j,t] = multi_student_t_rng(nu, mu[j,t],
					   multiply_lower_tri_self_transpose(L_H[j,t]))';
	log_lik[j,t] = multi_student_t_lpdf(rts[j,t] | nu, mu[j,t],
					    multiply_lower_tri_self_transpose(L_H[j,t]));
      }
    }
  }
  //
  for(j in 1:J){
    for(t in 1:T){
      //precision[j,t] = inverse(H[j,t]);
      H[j, t] = multiply_lower_tri_self_transpose(L_H[j,t]);
      pcor[j, t] = - cov2cor( inverse(multiply_lower_tri_self_transpose(L_H[j,t])) );
    }
  }
  Sfixed = invvec_to_corr( tanh(S_vec_fixed), nt);
  Sfixed2= invvec_to_corr( tanh(S_vec_fixed2), nt);
}

