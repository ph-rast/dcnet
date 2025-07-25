// Fixed effect vector
vector[nt] phi0_fixed;
vector[nt] phi0_fixed2;  // In case of S_pred != NULL 

cholesky_factor_corr[nt] phi0_L;
vector<lower=0>[nt] phi0_tau; // ranef SD for phi0
array[J] vector[nt] phi0_stdnorm; // Used to multiply with phi0_sd to obtain ranef on phi0

vector[nt*nt] vec_phi_fixed;

cholesky_factor_corr[nt*nt] phi_L;
vector<lower=0>[nt*nt] phi_tau; // ranef SD for phi0
array[J] vector[nt*nt] phi_stdnorm; // Used to multiply with phi0_sd to obtain ranef on phi0

// this creates uncorrelated random effects that are directly estimated
//vector[nt*nt] vec_phi_random[J];
//matrix<lower = -1 upper = 1>[nt,nt] theta;

