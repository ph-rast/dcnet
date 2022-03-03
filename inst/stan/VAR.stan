// VAR only model w. constant variance

data {
#include /data/data.stan
}

transformed data {                                  
  vector[nt] rts_m[J];
  vector[nt] rts_sd[J];
  int<lower=nt> Sdim = (nt + nt*nt) / 2 - 1; // Dimension of vec(S).
 
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

  // Residual covmat
  cov_matrix[nt] rescov;
  real< lower = 2 > nu; // nu for student_t
}

transformed parameters {
  // transform vec_phi to nt*nt parameter matrix
  //matrix<lower = -1, upper = 1>[nt,nt] phi;
  vector[nt] phi0[J]; // vector with fixed + random for intercept
  vector[nt*nt] phi[J]; // vectorized VAR parameter matrix, fixed + random

  vector[nt] mu[T,J];
  
  // Initialize t=1 across all subjects
  for( j in 1:J){
    mu[1,j,] = phi0_fixed;
  }

  // Beter to loop first through array, then matrix column and row 
  // iterations geq 2
  for (t in 2:T){
    for( j in 1:J){
      // Meanstructure model: DROP, if only VAR is allowed
#include /model_components/mu.stan
    }
  }
}

model {
  // priors
  // VAR
  phi0_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects
  phi_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects
  for(j in 1:J){
    phi0_stdnorm[J] ~ std_normal();
    phi_stdnorm[J] ~ std_normal();
  }
 
  // Prior for initial state
  //for(j in 1:J){
    rescov ~ wishart(nt + 1.0, diag_matrix(rep_vector(1.0, nt)) );
  //}
  // Prior on nu for student_t
  nu ~ normal( nt, 50 );
  phi0_fixed ~ multi_normal(rts_m, diag_matrix( rep_vector(1.0, nt) ) );
  vec_phi_fixed ~ normal(0, 0.5);
  
  // likelihood
  for( j in 1:J) {
    if ( distribution == 0 ) {
      for(t in 1:T){
	rts[t,j,] ~ multi_normal(mu[t,j,], rescov);
      }
    } else if ( distribution == 1 ) {
      for(t in 1:T){
	rts[t,j,] ~ multi_student_t(nu, mu[t,j,], rescov);
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
/* mar 2 */
}
