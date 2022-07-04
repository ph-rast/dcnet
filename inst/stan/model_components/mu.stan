/* if( meanstructure == 0 ){ */
/*   // Constant */
/*   mu[j,t ] = phi0_fixed; */
/*  } else if( meanstructure == 1 ){ */
/*   // arma11 */
/*   mu[j,t ] = phi0_fixed; //+ phi * rts[t-1, ] + theta * (rts[t-1, ] - mu[t-1,]) ; */
/*  } else if( meanstructure == 2 ){ */
  // VAR1

//phi0[j] = phi0_fixed + (diag_pre_multiply(phi0_tau, phi0_L)*phi0_stdnorm[j]);
//phi[j] = vec_phi_fixed + (diag_pre_multiply(phi_tau, phi_L)*phi_stdnorm[j]) ;

phi0[j] = phi0_fixed +
  (diag_pre_multiply(phi0_tau, diag_matrix(rep_vector(1.0, nt)))*phi0_stdnorm[j]);
phi[j] = vec_phi_fixed +
  (diag_pre_multiply(phi_tau, diag_matrix(rep_vector(1.0, nt*nt)))*phi_stdnorm[j]) ;


mu[j,t] = phi0[j] + to_matrix( phi[j], nt, nt ) * (rts[j,t-1]'- phi0[j]);//_fixed) ; //to_matrix fill by column major order
