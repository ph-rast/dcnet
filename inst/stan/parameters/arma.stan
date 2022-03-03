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

// this creates uncorrelated random effects that are directly estimated
//vector[nt*nt] vec_phi_random[J];
//matrix<lower = -1 upper = 1>[nt,nt] theta;

