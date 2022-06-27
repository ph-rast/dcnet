// DCC-Parameterization with random D components, random S fixed q's
functions {
#include /functions/jacobian.stan
#include /functions/invvec.stan
}

data {
#include /data/data.stan
}

transformed data {
  array[J] vector[nt] rts_m;
  array[J] vector[nt] rts_sd;
  int<lower=nt> Sdim = (nt + nt*nt) %/% 2 ; // Dimension of vec(S).
  
  // S = S_L*S_L | S_L is a cholesky factor
  // S_L = invec(S_Lv) | S_Lv is vec(S_L)
  // Note that S_Lv[1,1] is always 1, which leaves (nt + nt^2) /2 -1 parameters to be estimated
 
#include /transformed_data/xh_marker.stan
 
  if( meanstructure == 0 ){
    for ( j in 1:J ){
      for ( i in 1:nt ){
	rts_m[j] = rep_vector(mean(rts[j,i]),nt);
	rts_sd[j] =  rep_vector(sd(rts[j,i]),nt);
      }
    }
  } else if (meanstructure == 1 || meanstructure == 2 ){
    // set rts_m to first element in ts
    for ( j in 1:J ){
      rts_m[j] = rts[j,1]';
      rts_sd[j] = rep_vector(sd(rts[1,1]),nt);
    }
  }
}
parameters {
  // VAR parameters
#include /parameters/arma.stan
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
  
  /* // vector<lower = 0,  upper = 1 >[nt] a_h[Q]; */
  /* array[J,nt] simplex[Q] a_h_simplex; */
  /* //vector<lower=0, upper = 1>[nt] a_h_sum[J]; */
  /* array[J,nt] simplex[P] b_h_simplex; // Simplex for b_h within each timeseries */
  /* array[J] vector[nt] b_h_sum_s; */

    // GARCH q parameters
  //  real<lower=0, upper = 1 > a_q; // Define on log scale so that it can go below 0
  real l_a_q; // Define on log scale so that it can go below 0
  array[J] real l_a_q_r;
  real<lower=0> l_a_q_sigma; //random effect variance
  
  //real<lower=0, upper = (1 - a_q) > b_q; //
  real l_b_q; //
  array[J] real l_b_q_r;
  real<lower=0> l_b_q_sigma; //random effect variance
  
  corr_matrix[nt] S;  

  // R1 init
  array[J] corr_matrix[nt] R1_init;
  
  // Qr1 init
  array[J] cov_matrix[nt] Qr1_init;
  // D1 init
  array[J] vector<lower = 0>[nt] D1_init;
  // u1 init
  vector[nt] u1_init;

  real< lower = 2 > nu; // nu for student_t
}

transformed parameters {
  // transform vec_phi to nt*nt parameter matrix
  //matrix<lower = -1, upper = 1>[nt,nt] phi;
  array[J] vector[nt] phi0; // vector with fixed + random for intercept
  array[J] vector[nt*nt] phi; // vectorized VAR parameter matrix, fixed + random

  //Scale
  array[J] real<lower = 0, upper = 1> a_q;
  array[J] real<lower = 0, upper = 1> b_q;

  array[J] vector[nt] c_h;
  array[J] vector[nt] c_h_random; // variance on log metric
  //vector[nt] a_h[J];
  array[J] vector[nt] a_h_random; // variance on log metric
  array[J] vector[nt] b_h_random;

  array[J] vector<lower=0, upper = 1>[nt] a_h_sum;
  array[J] vector<lower=0, upper = 1>[nt] b_h_sum;
  
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
  
 
  // VAR phi parameter

  // Initialize t=1
  for( j in 1:J){
    
    mu[j,1] = rts[j, 1]'; //phi0_fixed;
    D[j,1] = sqrt( diagonal((rts[j]' * rts[j]) / (nt - 1)) );//D1_init[j];//
    u[j,1] = ( rts[j,1]' - mu[j,1] ) ./ D[j,1] ;//u1_init;
    Qr[j,1] = Qr1_init[j];//((rts[j]' * rts[j]) / (nt - 1));//
    Qr_sdi[j,1] = 1 ./ sqrt(diagonal(Qr[j,1])); //
    
    //R[j,1] = R1_init[j];//cov2cor((rts[j]' * rts[j]) / (nt - 1));
    // R[j,1]=diag_matrix(rep_vector(1.0, nt));
    R[j,1] = quad_form_diag(Qr[j,1], Qr_sdi[j,1]); //
    H[j,1] = quad_form_diag(R[j,1],  D[j,1]);

    ////////////
    // Ranefs //
    ////////////
    
    // R part
    a_q[j] =          1 ./ ( 1 + exp(-(l_a_q + l_a_q_r[j])) );
    b_q[j] = (1-a_q[j]) ./ ( 1 + exp(-(l_b_q + l_b_q_r[j])) );

    // D part
    // Should random effect corrleations be estimated:
    // Full estimation of corrmat with simplify_ch == 0
    // Estimate only diagonal with simplify_ch == 1
    if(simplify_ch==0){
      c_h_random[j] = (diag_pre_multiply(c_h_tau, c_h_L)*c_h_stdnorm[j]);
    } else if(simplify_ch == 1) {
      c_h_random[j] = (diag_pre_multiply(c_h_tau, diag_matrix(rep_vector(1.0, nt)))*c_h_stdnorm[j]);
    }
    c_h[j] = c_h_fixed + c_h_random[j];

    // Same approach here with simplify_ah
    if(simplify_ah==0){
      a_h_random[j] = diag_pre_multiply(a_h_tau, a_h_L)*a_h_stdnorm[j];
    } else if(simplify_ah == 1) {
      a_h_random[j] = diag_pre_multiply(a_h_tau, diag_matrix(rep_vector(1.0, nt)))*a_h_stdnorm[j];
    }
    
    // Bound sum of fixed and ranef between 0 and 1 with logistic function
    a_h_sum[j] = 1 ./ (1 + exp(-( a_h_fixed + a_h_random[j] )) );

    
    if(simplify_bh==0){      
      b_h_random[j] = (diag_pre_multiply(b_h_tau, b_h_L)*b_h_stdnorm[j]);
    } else if(simplify_bh == 1) {
      b_h_random[j] = (diag_pre_multiply(b_h_tau, diag_matrix(rep_vector(1.0, nt)))*b_h_stdnorm[j]);
    }
    // Bound sum of fixed and ranef between 0 and 1
    b_h_sum[j] = (1 - a_h_sum[j]) ./ (1 + exp(-(b_h_fixed + b_h_random[j])) );
    
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
      Qr[j,t ] = (1 - a_q[j] - b_q[j]) * S + a_q[j] * (u[j, t-1 ] * u[j, t-1 ]') + b_q[j] * Qr[j, t-1]; // S and UU' define dimension of Qr
      Qr_sdi[j, t] = 1 ./ sqrt(diagonal(Qr[j, t])) ; // inverse of diagonal matrix of sd's of Qr
      //    R[t,] = quad_form_diag(Qr[t,], inv(sqrt(diagonal(Qr[t,]))) ); // Qr_sdi[t,] * Qr[t,] * Qr_sdi[t,];
      R[j,t] = quad_form_diag(Qr[j,t], Qr_sdi[j,t]); // 
      H[j,t] = quad_form_diag(R[j,t],  D[j,t]);  // H = DRD;
    }
  }
}

model {
  // print("Upper Limits:", UPs);
  // priors
  l_a_q ~ normal(-1, 1);
  l_b_q ~ normal(1.7, 1);
  l_a_q_sigma ~ cauchy(0, 0.1);
  to_vector(l_a_q_r) ~ normal(0, l_a_q_sigma);
  l_b_q_sigma ~ cauchy(0, 0.1);
  to_vector(l_b_q_r) ~ normal(0, l_b_q_sigma);

  // VAR
  phi0_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects
  phi_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects
  //D part in DRD
  c_h_L ~ lkj_corr_cholesky(1);
  a_h_L ~ lkj_corr_cholesky(1);
  b_h_L ~ lkj_corr_cholesky(1);
  // R part in DRD

  phi0_tau ~ cauchy(0, 1); // SD for multiplication with cholesky phi0_L
  phi_tau ~ cauchy(0, 1); // SD for multiplication with cholesky phi0_L
  c_h_tau ~ cauchy(0, .5); // SD for c_h ranefs
  a_h_tau ~ cauchy(0, .5); // SD for c_h ranefs
  b_h_tau ~ cauchy(0, .5);

  
  // C
  to_vector(beta) ~ std_normal();
  to_vector(c_h_fixed) ~ std_normal();
  to_vector(a_h_fixed) ~ std_normal();
  to_vector(b_h_fixed) ~ std_normal();
  // Prior for initial state
  
  to_vector(u1_init) ~ std_normal();
  // Prior on nu for student_t
  //if ( distribution == 1 )
  nu ~ normal( nt, 50 );
  //to_vector(theta) ~ std_normal();
  //to_vector(phi) ~ std_normal();
  phi0_fixed ~ multi_normal(rts_m, diag_matrix( rep_vector(1.0, nt) ) );
  vec_phi_fixed ~ normal(0, 5);
  //  to_vector(a_h) ~ normal(0, .5);
  //to_vector(b_h) ~ normal(0, .5);
  S ~ lkj_corr( 1 );

  // likelihood
  for( j in 1:J) {
    R1_init[j] ~ lkj_corr( 1 );
    to_vector(D1_init[j]) ~ lognormal(-1, 1);
    Qr1_init[j] ~ wishart(nt + 1.0, diag_matrix(rep_vector(1.0, nt)) );
    // UL transform jacobian
    /* for(k in 1:nt) { */
    /*   ULs[j,k] ~ uniform(0, UPs[j,k]); // Truncation not needed. */
    /*   target += a_b_scale_jacobian(0.0, ULs[j,k], b_h_sum_s[j,k]); */
    /* } */
    phi0_stdnorm[J] ~ std_normal();
    phi_stdnorm[J] ~ std_normal();
    c_h_stdnorm[J] ~ std_normal();
    a_h_stdnorm[J] ~ std_normal();
    b_h_stdnorm[J] ~ std_normal();
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
/*   matrix[nt,T] rts_out[J]; */
/*   real log_lik[T]; */
/*   corr_matrix[nt] corH[T]; */
/*   // for the no-predictor case */
/*   vector<lower=0>[nt] c_h_var = exp(c_h); */
/*   // retrodict */
/* #include /generated/retrodict_H.stan */
/* 15:13 */
}

