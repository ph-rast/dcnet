int<lower=0> nobs;                // num of observations
int<lower=1> J;                   // number of groups or subjects
array[nobs] int<lower=1,upper=J> group; // vector with group ID
int<lower=2> T; // lengh of time series
int<lower=2> nt;    // number of time series
int<lower=1> Q; // MA component in MGARCH(P,Q), matrix A
int<lower=1> P; // AR component in MGARCH(P,Q), matrix B
array[J] matrix[T,nt] rts;  // multivariate time-series; array of length J of Txnt matrix
array[nobs] vector[nt] xC;  // TODO - match to rts if to be included // time-varying predictor for constant variance
int<lower=0, upper=1> distribution; // 0 = Normal; 1 = student_t
int<lower=0, upper=2> meanstructure; // Select model for location
