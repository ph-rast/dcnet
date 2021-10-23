// DCC-Parameterization
functions { 
#include /functions/cov2cor.stan
#include /functions/jacobian.stan
#include /functions/invvec.stan
}

data {
#include /data/data.stan
}

transformed data {                                  
  vector[nt] rts_m[J];
  vector[nt] rts_sd[J];
  int<lower=nt> Sdim = (nt + nt*nt) / 2 - 1; // Dimension of vec(S).
  // S = S_L*S_L | S_L is a cholesky factor
  // S_L = invec(S_Lv) | S_Lv is vec(S_L)
  // Note that S_Lv[1,1] is always 1, which leaves (nt - nt^2) /2 -1 parameters to be estimated
 
#include /transformed_data/xh_marker.stan
 
  if( meanstructure == 0 ){
    for ( j in 1:J ){
      for ( i in 1:nt ){
	rts_m[j,i] = mean(rts[,j,i]);
	rts_sd[j,i] =  sd(rts[,j,i]);
      }
    }
  } else if (meanstructure == 1 || meanstructure == 2 ){
      // set rts_m to first element in ts
    for ( j in 1:J ){
      for ( i in 1:nt ){
	rts_m[j,i] = rts[1,j,i];
	rts_sd[j,i] = sd(rts[,j,i]);
      }
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

  // Random effects fo c_h
  cholesky_factor_corr[nt] c_h_L;
  vector<lower=0>[nt] c_h_tau;
  vector[nt] c_h_stdnorm[J];
  // Random effects fo a_h
  cholesky_factor_corr[nt] a_h_L;
  vector<lower=0>[nt] a_h_tau; // ranef SD for phi0
  vector[nt] a_h_stdnorm[J]; // Used to multiply with phi0_sd to obtain ranef on phi0
  // Random effects fo b_h
  cholesky_factor_corr[nt] b_h_L;
  vector<lower=0>[nt] b_h_tau; // ranef SD for phi0
  vector[nt] b_h_stdnorm[J]; // Used to multiply with phi0_sd to obtain ranef on phi0


  
  // vector<lower = 0,  upper = 1 >[nt] a_h[Q];
  simplex[Q] a_h_simplex[J,nt];
  //vector<lower=0, upper = 1>[nt] a_h_sum[J];
  simplex[P] b_h_simplex[J,nt]; // Simplex for b_h within each timeseries
  vector[nt] b_h_sum_s[J];
  //vector[nt] b_h_sum_s[J]; // Unconstrained b_h_sum values. b_h[i] = U[i] b_h_simplex[i]; U[i] ~ U(0, 1 - sum(a_h[i]))
  // vector<lower = 0,  upper = 1 >[nt] b_h[P]; // TODO actually: 1 - a_h, across all Q and P...
  // GARCH q parameters 
  real<lower=0, upper = 1 > a_q; // 
  real<lower=0, upper = (1 - a_q) > b_q; //
  // Elements for random effects on S
  // Ranefs on S are put on VECH(S) lower tri, with dimension (nt + nt^2)/2
  cholesky_factor_corr[Sdim] S_L_R;  // Cholesky of random efx corrmat
  vector<lower=0>[Sdim] S_L_tau; //SD's for random efx
  vector[Sdim] S_L_stdnorm[J]; 
  vector[Sdim] S_Lv_fixed; // Vectorized fixed effect for S_L
  // Qr1 init
  cov_matrix[nt] Qr1_init;
  // D1 init
  vector<lower = 0>[nt] D1_init;
  // u1 init
  vector[nt] u1_init;

  real< lower = 2 > nu; // nu for student_t
}

transformed parameters {
  // transform vec_phi to nt*nt parameter matrix
  //matrix<lower = -1, upper = 1>[nt,nt] phi;
  vector[nt] phi0[J]; // vector with fixed + random for intercept
  vector[nt*nt] phi[J]; // vectorized VAR parameter matrix, fixed + random

  //Scale
  vector[nt] c_h[J]; 
  vector[nt] c_h_random[J]; // variance on log metric
  //vector[nt] a_h[J];
  vector[nt] a_h_random[J]; // variance on log metric
  vector[nt] b_h_random[J];

  vector<lower=0, upper = 1>[nt] a_h_sum[J];
  vector<lower=0, upper = 1>[nt] b_h_sum[J];
  
  // ranef vector for vec(S) part
  vector[Sdim] S_Lv_r[J];
  // S_Lv[j] = S_Lv_fixed + S_Lv_r[j]
  vector[Sdim] S_Lv[J];
  // S_Lv contains (nt+nt^2)/2-1 elements to be inv-vech into lower tri cholesky factor S_L
  // Note that S_L[1,1] = 1 as S = S_L*S_L is a corrmat
  //cholesky_factor_corr[nt] S_L[J];
  corr_matrix[nt] S[J];
  
  cov_matrix[nt] H[T,J];
  corr_matrix[nt] R[T,J];
  vector[nt] rr[T-1,J];
  vector[nt] mu[T,J];
  vector[nt] D[T,J];
  cov_matrix[nt] Qr[T,J];
  vector[nt] Qr_sdi[T,J];
  vector[nt] u[T,J];
  vector<lower = 0>[nt] vd[J];
  vector<lower = 0>[nt] ma_d[J];
  vector<lower = 0>[nt] ar_d[J];
  
  vector<lower=0, upper = 1>[nt] a_h[Q,J];// = simplex_to_bh(a_h_simplex[J], a_h_sum[J]);
  vector[nt] UPs[J];// = upper_limits(a_h[J]);
  vector[nt] ULs[J];// = raw_sum_to_b_h_sum(b_h_sum_s[J], UPs[J]);
  vector<lower = 0, upper = 1>[nt] b_h[P,J];// = simplex_to_bh(b_h_simplex[J], ULs[J]);

  // This vector is multiplied be simplex; has to be 0<x<1
  // take tanh( fixed + ranef )
  
  // VAR phi parameter
  //phi = 
  
  // Initialize t=1
  for( j in 1:J){
    for( q in 1:Q){
      a_h[q,j] = rep_vector(.5, nt);//simplex_to_bh(a_h_simplex[j], a_h_sum[j]);
    }
    UPs[j] = upper_limits(a_h[,j]);
    ULs[j] = raw_sum_to_b_h_sum(b_h_sum_s[j], UPs[j]);
    for( p in 1:P ){
      b_h[p,j] = rep_vector(.5, nt); //simplex_to_bh(b_h_simplex[j], ULs[j]);
    }
    
    mu[1,j,] = phi0_fixed;
    u[1,j] = u1_init;
    D[1,j] = D1_init;
    Qr[1,j] = Qr1_init;
    H[1,j] = Qr[1,j];
    R[1,j] = diag_matrix(rep_vector(1.0, nt));
    Qr_sdi[1,j] = rep_vector(1.0, nt);
  }
  
  // iterations geq 2
  for (t in 2:T){
    for( j in 1:J){
      // Meanstructure model: DROP, if only VAR is allowed
#include /model_components/mu.stan
      c_h_random[j] = (diag_pre_multiply(c_h_tau, c_h_L)*c_h_stdnorm[j]);
      c_h[j] = c_h_fixed + c_h_random[j];

      a_h_random[j] = (diag_pre_multiply(a_h_tau, a_h_L)*a_h_stdnorm[j]);
      // Bound sum of fixed and ranef between 0 and 1 with logistic function
      a_h_sum[j] = rep_vector(1.0, nt) ./ (1 + exp(-(a_h_fixed + a_h_random[j])) );

      b_h_random[j] = (diag_pre_multiply(b_h_tau, b_h_L)*b_h_stdnorm[j]);
      // Bound sum of fixed and ranef between 0 and 1
      b_h_sum[j] = rep_vector(1.0, nt) ./ (1 + exp(-(b_h_fixed + b_h_random[j])) );

      
      for(d in 1:nt){
	vd[j,d]   = 0.0;
	ma_d[j,d] = 0.0;
	ar_d[j,d]   = 0.0;
	// GARCH MA component
	for (q in 1:min( t-1, Q) ) {
	  rr[t-q, j, d] = square( rts[t-q, j, d] - mu[t-q, j, d] );
	  ma_d[j,d] = ma_d[j,d] + a_h[q, j, d]*rr[t-q, j, d] ;
	}
	// print("ma_d: ", "TS:", d, " Value:", ma_d[d], " T:", t);
	// GARCH AR component
	for (p in 1:min( t-1, P) ) {
	  ar_d[j,d] = ar_d[j,d] + b_h[p, j, d]*D[t-p, j, d]^2;
	}
	// print("ar_d: ", "TS:", d, " Value:", ar_d[d], " T:", t);
	// Predictor on diag (given in xC)
	if ( xC_marker >= 1) {
	  vd[j,d] = exp( c_h[j,d] + beta[d] * xC[t, d] ) + ma_d[j,d] + ar_d[j,d];
	} else if ( xC_marker == 0) {
	  vd[j,d] = exp( c_h[j,d] )  + ma_d[j,d] + ar_d[j,d];
	}

	D[t, j, d] = sqrt( vd[j,d] );
      }
      u[t,j,] = diag_matrix(D[t,j,]) \ (rts[t,j,]'- mu[t,j,]); // cf. comment about taking inverses in stan manual p. 482 re:Inverses - inv(D)*y = D \ a

      // All output is in vectorized form
      S_Lv_r[j] = (diag_pre_multiply(S_L_tau, S_L_R)*S_L_stdnorm[j]);
      S_Lv[j] = S_Lv_fixed + S_Lv_r[j]; // keep within -1; 1
      // S_Lv is vectorized - invvec now:
      S[j] = invvec_to_corr(S_Lv[j], nt);
      
      Qr[t,j,] = (1 - a_q - b_q) * S[j] + a_q * (u[t-1,j,] * u[t-1,j,]') + b_q * Qr[t-1,j,]; // S and UU' define dimension of Qr
      Qr_sdi[t,j,] = 1 ./ sqrt(diagonal(Qr[t,j,])); // inverse of diagonal matrix of sd's of Qr
      //    R[t,] = quad_form_diag(Qr[t,], inv(sqrt(diagonal(Qr[t,]))) ); // Qr_sdi[t,] * Qr[t,] * Qr_sdi[t,];
      R[t,j] = quad_form_diag(Qr[t,j,], Qr_sdi[t,j,]); // Qr_sdi[t,] * Qr[t,] * Qr_sdi[t,];
      H[t,j] = quad_form_diag(R[t,j],     D[t,j]);  // H = DRD; 
    }
  }
}
model {
  // print("Upper Limits:", UPs);
  // UL transform jacobian
  for(k in 1:nt) {
    for(j in 1:J) {
      ULs[j,k] ~ uniform(0, UPs[j,k]); // Truncation not needed.
      target += a_b_scale_jacobian(0.0, ULs[j,k], b_h_sum_s[j,k]);
    }
  }

  // priors
  // VAR
  phi0_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects
  phi_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects
  //D part in DRD
  c_h_L ~ lkj_corr_cholesky(1);
  a_h_L ~ lkj_corr_cholesky(1);
  b_h_L ~ lkj_corr_cholesky(1);
  // R part in DRD
  S_L_R ~ lkj_corr_cholesky(1);
  phi0_tau ~ cauchy(0, 2); // SD for multiplication with cholesky phi0_L
  phi_tau ~ cauchy(0, 2); // SD for multiplication with cholesky phi0_L
  c_h_tau ~ cauchy(0, 2); // SD for c_h ranefs
  a_h_tau ~ cauchy(0, 2); // SD for c_h ranefs
  b_h_tau ~ cauchy(0, 2);
  S_L_tau ~ cauchy(0, 2); 
  for(j in 1:J){
    phi0_stdnorm[J] ~ std_normal();
    phi_stdnorm[J] ~ std_normal();
    c_h_stdnorm[J] ~ std_normal();
    a_h_stdnorm[J] ~ std_normal();
    b_h_stdnorm[J] ~ std_normal();
    S_L_stdnorm[J] ~ std_normal();
  }


  
  // C
  to_vector(beta) ~ std_normal();
  to_vector(c_h_fixed) ~ std_normal();
  to_vector(a_h_fixed) ~ std_normal();
  to_vector(b_h_fixed) ~ std_normal();
  // Prior for initial state
  Qr1_init ~ wishart(nt + 1.0, diag_matrix(rep_vector(1.0, nt)) );
  to_vector(D1_init) ~ lognormal(-1, 1);
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
  //  S ~ lkj_corr( 1 );

  // likelihood
  for( j in 1:J) {
    if ( distribution == 0 ) {
      for(t in 1:T){
	rts[t,j,] ~ multi_normal(mu[t,j,], H[t,j,]);
      }
    } else if ( distribution == 1 ) {
      for(t in 1:T){
	rts[t,j,] ~ multi_student_t(nu, mu[t,j,], H[t,j,]);
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
