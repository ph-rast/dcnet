// VAR only model w. constant variance

data {
#include /data/data.stan
}

transformed data {                                  
  array[J] vector[nt] rts_m;
  array[J] vector[nt] rts_sd;
  int<lower=nt> Sdim = (nt + nt*nt) %/% 2 ; // Dimension of vec(S).
 
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
  
  // Residual covmat
  cov_matrix[nt] rescov;
  real< lower = 2 > nu; // nu for student_t
}

transformed parameters {
  // transform vec_phi to nt*nt parameter matrix
  //matrix<lower = -1, upper = 1>[nt,nt] phi;
  array[J] vector[nt] phi0; // vector with fixed + random for intercept
  array[J] vector[nt*nt] phi; // vectorized VAR parameter matrix, fixed + random

  array[J,T] vector[nt] mu;
  
  // Initialize t=1 across all subjects
  for( j in 1:J){
    mu[j,1] = rts[j, 1]';
  }

  // Beter to loop first through array, then matrix column and row 
  // iterations geq 2
  for( j in 1:J ){
    for( t in 2:T ){
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
    phi0_stdnorm[j] ~ std_normal();
    phi_stdnorm[j] ~ std_normal();
  }
  // Prior for initial state
  rescov ~ wishart(nt + 1.0, diag_matrix(rep_vector(1.0, nt)) );
  // Prior on nu for student_t
  nu ~ normal( nt, 50 );
  //  phi0_fixed ~ multi_normal(rts_m[1,1], diag_matrix( rep_vector(1.0, nt) ) ); //REVISIT rts[1,1]
  phi0_fixed ~ normal(0, 10 );
  vec_phi_fixed ~ normal(0, 0.5);
  
  // likelihood 
  for( j in 1:J) { 
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
