if (S_pred[j,t] == 0){
  phi0[j] = phi0_fixed + phi0_tau .* phi0_stdnorm[j];
 } else if (S_pred[j,t] == 1){
  phi0[j] = phi0_fixed + phi0_fixed2 + 
    (diag_pre_multiply(phi0_tau, diag_matrix(rep_vector(1.0, nt)))*phi0_stdnorm[j]);
 }
  
phi[j] = vec_phi_fixed + phi_tau .* phi_stdnorm[j] ;

mu[j,t] = phi0[j] + to_matrix( phi[j], nt, nt ) * (rts[j,t-1]'- phi0[j]); //to_matrix fill by column major order
