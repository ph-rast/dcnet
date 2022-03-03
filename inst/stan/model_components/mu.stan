if( meanstructure == 0 ){
  // Constant
  mu[t,j, ] = phi0_fixed;
 } else if( meanstructure == 1 ){
  // arma11
  mu[t,j, ] = phi0_fixed; //+ phi * rts[t-1, ] + theta * (rts[t-1, ] - mu[t-1,]) ;
 } else if( meanstructure == 2 ){
  // VAR1
  phi0[j] = phi0_fixed + (diag_pre_multiply(phi0_tau, phi0_L)*phi0_stdnorm[j]);
  phi[j] = vec_phi_fixed + (diag_pre_multiply(phi_tau, phi_L)*phi_stdnorm[j]) ;
  mu[t,j, ] = phi0[j] + to_matrix( phi[j], nt, nt ) * (rts[,t-1,j]'- phi0_fixed); //phi0[j]);
    //  to_matrix( (vec_phi_fixed + vec_phi_random[j,]), nt, nt ) * (rts[t-1,j,]'-phi0[j]);
  //phi * rts[t-1,j,]'; // vectorize phi, add ranefs and put it back together
 }

// could use tanh(fixed + random) to keep it all within -1;1 

