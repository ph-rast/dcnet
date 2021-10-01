// DCC-Parameterization
functions { 
#include /functions/cov2cor.stan
#include /functions/jacobian.stan
}

data {
#include /data/data.stan
}

transformed data {                                  
  vector[nt] rts_m[J];
  vector[nt] rts_sd[J];
 
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
  vector[nt] c_h; // variance on log metric 
  // vector<lower = 0,  upper = 1 >[nt] a_h[Q];
  simplex[Q] a_h_simplex[nt];
  vector<lower=0, upper = 1>[nt] a_h_sum;
  simplex[P] b_h_simplex[nt]; // Simplex for b_h within each timeseries
  vector[nt] b_h_sum_s; // Unconstrained b_h_sum values. b_h[i] = U[i] b_h_simplex[i]; U[i] ~ U(0, 1 - sum(a_h[i]))
  // vector<lower = 0,  upper = 1 >[nt] b_h[P]; // TODO actually: 1 - a_h, across all Q and P...
  // GARCH q parameters 
  real<lower=0, upper = 1 > a_q; // 
  real<lower=0, upper = (1 - a_q) > b_q; //
  corr_matrix[nt] S;  // DCC keeps this constant 
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
  
  cov_matrix[nt] H[T];
  corr_matrix[nt] R[T];
  vector[nt] rr[T-1];
  vector[nt] mu[T,J];
  vector[nt] D[T];
  cov_matrix[nt] Qr[T];
  vector[nt] Qr_sdi[T];
  vector[nt] u[T];
  real<lower = 0> vd[nt];
  real<lower = 0> ma_d[nt];
  real<lower = 0> ar_d[nt];  
  vector<lower=0, upper = 1>[nt] a_h[Q] = simplex_to_bh(a_h_simplex, a_h_sum);
  vector[nt] UPs = upper_limits(a_h);
  vector[nt] ULs = raw_sum_to_b_h_sum(b_h_sum_s, UPs);
  vector<lower = 0, upper = 1>[nt] b_h[P] = simplex_to_bh(b_h_simplex, ULs);

  // VAR phi parameter
  //phi = 
  
  // Initialize t=1
  for( j in 1:J){
    mu[1,j,] = phi0_fixed;
  }
  u[1,] = u1_init;
  D[1,] = D1_init;
  Qr[1,] = Qr1_init;
  H[1] = Qr[1,];
  R[1] = diag_matrix(rep_vector(1.0, nt));
  Qr_sdi[1] = rep_vector(1.0, nt);

  // iterations geq 2
  for (t in 2:T){
    for( j in 1:J){
      // Meanstructure model: DROP, if only VAR is allowed
#include /model_components/mu.stan 

      for(d in 1:nt){
	vd[d] = 0.0;
	ma_d[d] = 0.0;
	ar_d[d] = 0.0;
	// GARCH MA component
	for (q in 1:min( t-1, Q) ) {
	  rr[t-q, d] = square( rts[t-q, j, d] - mu[t-q, j, d] );
	  ma_d[d] = ma_d[d] + a_h[q, d]*rr[t-q, d] ;
	}
	// print("ma_d: ", "TS:", d, " Value:", ma_d[d], " T:", t);
	// GARCH AR component
	for (p in 1:min( t-1, P) ) {
	  ar_d[d] = ar_d[d] + b_h[p, d]*D[t-p, d]^2;
	}
	// print("ar_d: ", "TS:", d, " Value:", ar_d[d], " T:", t);
	// Predictor on diag (given in xC)
	if ( xC_marker >= 1) {
	  vd[d] = exp( c_h[d] + beta[d] * xC[t, d] ) + ma_d[d] + ar_d[d];
	} else if ( xC_marker == 0) {
	  vd[d] = exp( c_h[d] )  + ma_d[d] + ar_d[d];
	}
	// print("c_h: ", "TS: ", d, " Value:", c_h[d]);

	D[t, d] = sqrt( vd[d] );
      }
      u[t,] = diag_matrix(D[t,]) \ (rts[t,j,]'- mu[t,j,]); // cf. comment about taking inverses in stan manual p. 482 re:Inverses - inv(D)*y = D \ a
      // print("u: ", "Value: ", u[t], "T:", t);
      Qr[t,] = (1 - a_q - b_q) * S + a_q * (u[t-1,] * u[t-1,]') + b_q * Qr[t-1,]; // S and UU' define dimension of Qr
      Qr_sdi[t,] = 1 ./ sqrt(diagonal(Qr[t,])); // inverse of diagonal matrix of sd's of Qr
      //    R[t,] = quad_form_diag(Qr[t,], inv(sqrt(diagonal(Qr[t,]))) ); // Qr_sdi[t,] * Qr[t,] * Qr_sdi[t,];
      R[t,] = quad_form_diag(Qr[t,], Qr_sdi[t,]); // Qr_sdi[t,] * Qr[t,] * Qr_sdi[t,];
      H[t,] = quad_form_diag(R[t,],     D[t,]);  // H = DRD; 
    }
  }
}
model {
  // print("Upper Limits:", UPs);
  // UL transform jacobian
  for(k in 1:nt) {
    ULs[k] ~ uniform(0, UPs[k]); // Truncation not needed.
    target += a_b_scale_jacobian(0, ULs[k], b_h_sum_s[k]);
  }

  // priors
  // VAR
  phi0_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects
  phi_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects
  phi0_tau ~ cauchy(0, 2); // SD for multiplication with cholesky phi0_L
  phi_tau ~ cauchy(0, 2); // SD for multiplication with cholesky phi0_L
  for(j in 1:J){
    phi0_stdnorm[J] ~ std_normal();
    phi_stdnorm[J] ~ std_normal();
  }
  // C
  to_vector(beta) ~ std_normal();
  to_vector(c_h) ~ std_normal();
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
  S ~ lkj_corr( 1 );

  // likelihood
  for( j in 1:J) {
    if ( distribution == 0 ) {
      for(t in 1:T){
	rts[t,j,] ~ multi_normal(mu[t,j,], H[t,]);
      }
    } else if ( distribution == 1 ) {
      for(t in 1:T){
	rts[t,j,] ~ multi_student_t(nu, mu[t,j,], H[t,]);
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
/* 00:42 */
}
