int<lower=0> nobs;                // num of observations
int<lower=1> J;                   // number of groups or subjects
array[nobs] int<lower=1,upper=J> group; // vector with group ID
int<lower=2> T; // lengh of time series
int<lower=2> nt;    // number of time series
int<lower=1> Q; // MA component in MGARCH(P,Q), matrix A
int<lower=1> P; // AR component in MGARCH(P,Q), matrix B
real lbound; //lower bound of observed variables
real ubound; //upper bound of observed variables
array[J] matrix<lower=lbound,upper=ubound>[T,nt] rts;  // multivariate time-series; array of length J of Txnt matrix
array[nobs] vector[nt] xC;  // TODO - match to rts if to be included // time-varying predictor for constant variance
int<lower=0, upper=1> distribution; // 0 = Normal; 1 = student_t
int<lower=0, upper=2> meanstructure; // Select model for location
int<lower=0, upper=1> simplify_ch; // 0 = random effects correlations among nt's, 1 = only ranefs without corrs among nt's
int<lower=0, upper=1> simplify_ah; // 0 = random effects correlations among nt's, 1 = only ranefs without corrs among nt's
int<lower=0, upper=1> simplify_bh; // 0 = random effects correlations among nt's, 1 = only ranefs without corrs among nt's
array[J] vector<lower=0, upper=1>[T] S_pred; // predictor for S
vector[nt] grand_mean;
int<lower=1> grainsize;  // e.g., 1, 2, or larger depending on your hardware
