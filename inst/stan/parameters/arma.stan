vector[nt] phi0; // Fixed effect vector

cholesky_factor_corr[nt] phi0_L;
vector<lower=0>[nt] phi0_tau; // ranef SD for phi0
vector[nt] phi0_stdnorm[J]; // Used to multiply with phi0_sd to obtain ranef on phi0
  
// gqs() cant not deal with ? yet - as long as that's not fixed
// estimate phi and theta anyways
//matrix[meanstructure ? nt : 0, meanstructure ? nt : 0 ] phi;
//matrix[meanstructure ? nt : 0, meanstructure ? nt : 0 ] theta;
matrix<lower = -1, upper = 1>[nt,nt] phi;
matrix<lower = -1, upper = 1>[nt,nt] theta;
