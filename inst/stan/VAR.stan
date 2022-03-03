// VAR only model w. constant variance

data {
  //#include /data/data.stan
int<lower=0> nobs;                // num of observations
int<lower=1> J;                   // number of groups or subjects
int<lower=1,upper=J> group[nobs]; // vector with group ID
int<lower=2> T; // lengh of time series
int<lower=2> nt;    // number of time series
int<lower=1> Q; // MA component in MGARCH(P,Q), matrix A
int<lower=1> P; // AR component in MGARCH(P,Q), matrix B
matrix[50,3] rts[5];  // multivariate time-series; array of length J of Txnt matrix
vector[nt] xC[nobs];  // TODO - match to rts if to be included // time-varying predictor for constant variance
int<lower=0, upper=1> distribution; // 0 = Normal; 1 = student_t
int<lower=0, upper=2> meanstructure; // Select model for location

}

transformed data {                                  
  vector[nt] rts_m[J];
  vector[nt] rts_sd[J];
  int<lower=nt> Sdim = (nt + nt*nt) / 2 - 1; // Dimension of vec(S).
 
  //#include /transformed_data/xh_marker.stan
// Check whether xC contains a predictor or not.
matrix[nt, nt] xC_m[T];
int<lower = 0> xC_marker = 0;
real<lower = 0> cp;

for( t in 1:T ){
  xC_m[t] = diag_matrix( xC[t] );
  // add a variable that notes if xC is null or actually a predictor
  cp = sum( xC_m[t]' * xC_m[t] );
  if( cp != 0 )
    xC_marker = xC_marker + 1;
 }

  
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
  //#include /parameters/arma.stan
// Fixed effect vector
vector[nt] phi0_fixed; 
cholesky_factor_corr[nt] phi0_L;
vector<lower=0>[nt] phi0_tau; // ranef SD for phi0
vector[nt] phi0_stdnorm[J]; // Used to multiply with phi0_sd to obtain ranef on phi0
  
// gqs() cant not deal with ? yet - as long as that's not fixed
// estimate phi and theta anyways
//matrix[meanstructure ? nt : 0, meanstructure ? nt : 0 ] phi;
//matrix[meanstructure ? nt : 0, meanstructure ? nt : 0 ] theta;

//matrix<lower = -1, upper = 1>[nt,nt] phi; // TODO lower and upper may need to be removed when adding ranefs
// Create vecotr for phi here as parameter, then combine it to matrix in transformed parameters block
vector[nt*nt] vec_phi_fixed;
cholesky_factor_corr[nt*nt] phi_L;
vector<lower=0>[nt*nt] phi_tau; // ranef SD for phi0
vector[nt*nt] phi_stdnorm[J]; // Used to multiply with phi0_sd to obtain ranef on phi0
  

  // Residual covmat
  cov_matrix[nt] rescov;
  real< lower = 2 > nu; // nu for student_t
}

transformed parameters {
  // transform vec_phi to nt*nt parameter matrix
  //matrix<lower = -1, upper = 1>[nt,nt] phi;
  vector[nt] phi0[J]; // vector with fixed + random for intercept
  vector[nt*nt] phi[J]; // vectorized VAR parameter matrix, fixed + random

  matrix[T,nt] mu[J];
  
  // Initialize t=1 across all subjects
  for( j in 1:J){
    mu[j,1] = phi0_fixed';
  }

  // Beter to loop first through array, then matrix column and row 
  // iterations geq 2
  for( j in 1:J ){
    for( t in 2:T ){
      //#include /model_components/mu.stan
      phi0[j] = phi0_fixed + (diag_pre_multiply(phi0_tau, phi0_L)*phi0_stdnorm[j]);
phi[j] = vec_phi_fixed + (diag_pre_multiply(phi_tau, phi_L)*phi_stdnorm[j]) ;
mu[j,t] = ( phi0[j] + to_matrix( phi[j], nt, nt ) * (rts[j,t-1]'- phi0_fixed) )'; //phi0[j]);

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
/*   matrix[nt,T] rts_out[J]; */
/*   real log_lik[T]; */
/*   corr_matrix[nt] corH[T]; */
/*   // for the no-predictor case */
/*   vector<lower=0>[nt] c_h_var = exp(c_h); */
/*   // retrodict */
/* #include /generated/retrodict_H.stan */
}
