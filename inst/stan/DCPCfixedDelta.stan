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
  /* vector[nt] c_h_fixed; // variance on log metric */
  /* vector[nt] a_h_fixed; */
  /* vector[nt] b_h_fixed; */

  /* // Random effects for c_h */
  /* cholesky_factor_corr[nt] c_h_L; */
  /* vector<lower=0>[nt] c_h_tau; */
  /* array[J] vector[nt] c_h_stdnorm; */
  
  /* // Random effects fo a_h */
  /* cholesky_factor_corr[nt] a_h_L; */
  /* vector<lower=0>[nt] a_h_tau; // ranef SD for phi0 */
  /* array[J] vector[nt] a_h_stdnorm; // Used to multiply with phi0_sd to obtain ranef on phi0 */

  /* // Random effects fo b_h */
  /* cholesky_factor_corr[nt] b_h_L; */
  /* vector<lower=0>[nt] b_h_tau; // ranef SD for phi0 */
  /* array[J] vector[nt] b_h_stdnorm; // Used to multiply with phi0_sd to obtain ranef on phi0 */
  
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
  corr_matrix[nt] S2;

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

  /* array[J] vector[nt] c_h; */
  /* array[J] vector[nt] c_h_random; // variance on log metric */
  /* //vector[nt] a_h[J]; */
  /* array[J] vector[nt] a_h_random; // variance on log metric */
  /* array[J] vector[nt] b_h_random; */

  /* array[J] vector<lower=0, upper = 1>[nt] a_h_sum; */
  /* array[J] vector<lower=0, upper = 1>[nt] b_h_sum; */
  
  array[J,T] cov_matrix[nt] H;
  array[J,T] corr_matrix[nt] R;  // R is part cor
  array[J,T] corr_matrix[nt] IR; // the I-R marix in D (I-R)^-1 D
  array[J,T-1] vector[nt] rr;
  array[J,T] vector[nt] mu;
  vector<lower = 0>[nt] D;
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
    //D[j,1] = sqrt( diagonal((rts[j]' * rts[j]) / (nt - 1)) );//D1_init[j];//
    //u[j,1] = ( rts[j,1]' - mu[j,1] ) ./ D[j,1] ;//u1_init;
    u[j,1] = u1_init;    
    Qr[j,1] = Qr1_init[j];//

    Qr_sdi[j,1] = 1 ./ sqrt(diagonal(Qr[j,1])); //
    
    R[j,1] = diag_matrix(rep_vector(1.0, nt));//quad_form_diag(Qr[j,1], Qr_sdi[j,1]); //
    H[j,1] = diag_matrix(rep_vector(1.0, nt));//quad_form_diag(R[j,1],  D[j,1]);

    IR[j,1] = diag_matrix(rep_vector(1.0, nt)); // Fill first matrix

    ////////////
    // Ranefs //
    ////////////
    
    // R part
    a_q[j] =          1 ./ ( 1 + exp(-(l_a_q + l_a_q_r[j])) );
    b_q[j] = (1-a_q[j]) ./ ( 1 + exp(-(l_b_q + l_b_q_r[j])) );

    
    /////////////////////
    // Model Structure //
    /////////////////////
    
    for (t in 2:T){
      // Meanstructure model: DROP, if only VAR is allowed
#include /model_components/mu.stan     
      
      // Extract SD of H matrix as in H = SD*R*SD
      u[j,t,] = diag_matrix( 1 ./ sqrt(diagonal(H[j,t])) ) \ (rts[j,t]' - mu[j,t]);

      if (S_pred[j,t] == 0){
	Qr[j,t ] = (1 - a_q[j] - b_q[j]) * S + a_q[j] * (u[j, t-1 ] * u[j, t-1 ]') + b_q[j] * Qr[j, t-1]; // S and UU' define dimension of Qr
      } else if (S_pred[j,t] == 1){
	Qr[j,t ] = (1 - a_q[j] - b_q[j]) * S2 + a_q[j] * (u[j, t-1 ] * u[j, t-1 ]') + b_q[j] * Qr[j, t-1];
      }
      
      Qr_sdi[j, t] = 1 ./ sqrt(diagonal(Qr[j, t])) ; // inverse of diagonal matrix of sd's of Qr
      //    R[t,] = quad_form_diag(Qr[t,], inv(sqrt(diagonal(Qr[t,]))) ); // Qr_sdi[t,] * Qr[t,] * Qr_sdi[t,];
      R[j,t] = quad_form_diag(Qr[j,t], Qr_sdi[j,t]); //Partial Correlatin Matrix
      IR[j,t] = inverse( -R[j,t] + diag_matrix(rep_vector(2.0, nt)) ); //Instead of I - R, do -R + 2, so that diag of R = 1 instead of 0; This is a corrmat, I think...
      H[j,t] = quad_form_diag( IR[j,t],  D);  //  D here must contain the inverted sqrt'd precision elements; 1/sqrt(diag(theta))
      //H[j,t] = quad_form_diag(R[j,t],  D[j,t]);  // H = DRD;
    }
  }
}

model {
  // print("Upper Limits:", UPs);
  // priors
  l_a_q ~ normal(-1, 1);
  l_b_q ~ normal(2, 1);
  l_a_q_sigma ~ cauchy(0, 0.1);
  to_vector(l_a_q_r) ~ normal(0, l_a_q_sigma);
  l_b_q_sigma ~ cauchy(0, 0.1);
  to_vector(l_b_q_r) ~ normal(0, l_b_q_sigma);

  // VAR
  phi0_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects
  phi_L ~ lkj_corr_cholesky(1); // Cholesky of location random intercept effects

  // R part in DRD

  phi0_tau ~ cauchy(0, 1); // SD for multiplication with cholesky phi0_L
  phi_tau ~ cauchy(0, 1); // SD for multiplication with cholesky phi0_L
    
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
  S ~ lkj_corr( 10 );
  S2 ~ lkj_corr( 10 );

  to_vector(D) ~ lognormal(-1, 1);
  
  // likelihood
  for( j in 1:J) {
    //R1_init[j] ~ lkj_corr( 1 );
    Qr1_init[j] ~ wishart(nt + 1.0, diag_matrix(rep_vector(1.0, nt)) );
    phi0_stdnorm[J] ~ std_normal();
    phi_stdnorm[J] ~ std_normal();
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
#include /generated/retrodict.stan
}

