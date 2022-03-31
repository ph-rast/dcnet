// DCC-Parameterization
functions {
  //#include /functions/cov2cor.stan
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
 
  // GARCH q parameters
  //  real<lower=0, upper = 1 > a_q; // Define on log scale so that it can go below 0
  real l_a_q; // Define on log scale so that it can go below 0
  array[J] real l_a_q_r;
  //real<lower=0, upper = (1 - a_q) > b_q; //
  real l_b_q; //
  array[J] real l_b_q_r;
  // Elements for random effects on S
  // Ranefs on S are put on VECH(S) lower tri, with dimension (nt + nt^2)/2
  cholesky_factor_corr[Sdim] S_L_R;  // Cholesky of random efx corrmat
  vector<lower=0>[Sdim] S_L_tau; //SD's for random efx
  array[J] vector[Sdim] S_L_stdnorm;
  vector[Sdim] S_Lv_fixed; // Vectorized fixed effect for S_L
  // Qr1 init
  cov_matrix[nt] Qr1_init;
  // D1 init
  //vector<lower = 0>[nt] D1_init;
  // u1 init
  vector[nt] u1_init;

  real< lower = 2 > nu; // nu for student_t

  vector[nt] D;

}

transformed parameters {
  // transform vec_phi to nt*nt parameter matrix
  //matrix<lower = -1, upper = 1>[nt,nt] phi;
  array[J] vector[nt] phi0; // vector with fixed + random for intercept
  array[J] vector[nt*nt] phi; // vectorized VAR parameter matrix, fixed + random

  //Scale
  array[J] real<lower = 0, upper = 1> a_q;
  array[J] real<lower = 0, upper = 1> b_q;
  
  // ranef vector for vec(S) part
  array[J] vector[Sdim] S_Lv_r;
  // S_Lv[j] = S_Lv_fixed + S_Lv_r[j]
  array[J] vector[Sdim] S_Lv;
  // S_Lv contains (nt+nt^2)/2-1 elements to be inv-vech into lower tri cholesky factor S_L
  // Note that S_L[1,1] = 1 as S = S_L*S_L is a corrmat
  //cholesky_factor_corr[nt] S_L[J];
  array[J] corr_matrix[nt] S;
  
  array[J,T] cov_matrix[nt] H;
  array[J,T] corr_matrix[nt] R;
  array[J,T-1] vector[nt] rr;
  array[J,T] vector[nt] mu;
  array[J,T] cov_matrix[nt] Qr;
  array[J,T] vector[nt] Qr_sdi;
  array[J,T] vector[nt] u;
 
  // This vector is multiplied be simplex; has to be 0<x<1
  // take tanh( fixed + ranef )
  
  // VAR phi parameter
  //phi =
  
  // Initialize t=1
  for( j in 1:J){
    
    mu[j,1] = phi0_fixed;
    u[j,1] = u1_init;
    
    Qr[j,1] = Qr1_init;
    //H[j,1] = Qr[j,1];
    Qr_sdi[j,1] = 1 ./ sqrt(diagonal(Qr[j,1])); //
    R[j,1] = quad_form_diag(Qr[j,1], Qr_sdi[j,1]); //
    H[j,1] = quad_form_diag(R[j,1],  D );
    
    a_q[j] =          1 ./ ( 1 + exp(-(l_a_q + l_a_q_r[j])) );
    b_q[j] = (1-a_q[j]) ./ ( 1 + exp(-(l_b_q + l_b_q_r[j])) );
    
    for (t in 2:T){
      // Meanstructure model: DROP, if only VAR is allowed
#include /model_components/mu.stan
      
      u[j,t,] = diag_matrix(D) \ (rts[j,t]'- mu[j,t]) ; // cf. comment about taking inverses in stan manual p. 482 re:Inverses - inv(D)*y = D \ a

      // All output is in vectorized form
      S_Lv_r[j] = (diag_pre_multiply(S_L_tau, S_L_R)*S_L_stdnorm[j]); // SD metric
      // could take atanh here of fixed + random to create S_Lv and then have one element less
      // to estimate as cor[1,1] element in chol is always 1 for correlations
      S_Lv[j] = S_Lv_fixed + S_Lv_r[j]; //S_Lv(_fixed) is on cholesky L cov metric 
      // S_Lv is vectorized - invvec now and return cor: 
      S[j] = invvec_to_corr(S_Lv[j], nt);
      
      Qr[j,t ] = (1 - a_q[j] - b_q[j]) * S[j] + a_q[j] * (u[j, t-1 ] * u[j, t-1 ]') + b_q[j] * Qr[j, t-1]; // S and UU' define dimension of Qr
      Qr_sdi[j, t] = 1 ./ sqrt(diagonal(Qr[j, t])) ; // inverse of diagonal matrix of sd's of Qr
      //    R[t,] = quad_form_diag(Qr[t,], inv(sqrt(diagonal(Qr[t,]))) ); // Qr_sdi[t,] * Qr[t,] * Qr_sdi[t,];
      R[j,t] = quad_form_diag(Qr[j,t], Qr_sdi[j,t]); // 
      H[j,t] = quad_form_diag(R[j,t],  D );  // H = DRD;
    }
  }
}
model {

  // priors
  l_a_q ~ normal(-3, 1);
  l_b_q ~ normal(1, 1);
  to_vector(l_a_q_r) ~ std_normal();
  to_vector(l_b_q_r) ~ std_normal();
  // VAR
  phi0_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects
  phi_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects
  //D part in DRD
  // R part in DRD
  S_L_R ~ lkj_corr_cholesky(1);
  phi0_tau ~ cauchy(0, 1); // SD for multiplication with cholesky phi0_L
  phi_tau ~ cauchy(0, 1); // SD for multiplication with cholesky phi0_L
  S_L_tau ~ cauchy(0, .25);

  for(j in 1:J){
    phi0_stdnorm[J] ~ std_normal();
    phi_stdnorm[J] ~ std_normal();
    S_L_stdnorm[J] ~ std_normal();
  }
 
  // C
  // to_vector(beta) ~ std_normal();
  // Prior for initial state
  Qr1_init ~ wishart(nt + 1.0, diag_matrix(rep_vector(1.0, nt)) );
  to_vector(D) ~ cauchy(0, 1);
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

