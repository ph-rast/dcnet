if( meanstructure == 0 ){
  // Constant
  mu[t,j, ] = phi0;
 } else if( meanstructure == 1 ){
  // arma11
  mu[t,j, ] = phi0; //+ phi * rts[t-1, ] + theta * (rts[t-1, ] - mu[t-1,]) ;
 } else if( meanstructure == 2 ){
  // VAR1
  mu[t,j, ] = phi0 + (diag_pre_multiply(phi0_tau, phi0_L)*phi0_stdnorm[j]) + phi * rts[t-1,j,]';
 }
