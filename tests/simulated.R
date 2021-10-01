library(dcnet )

## Create data:
N <- 20
tl <- 20
nts <- 3
simdat <- dcnet:::.simuDCC(tslength = tl,  N = N,  n_ts = nts,  ranef_sd_S = 0.0001 )
simdat[[1]]

## Generated TS for person 1: simdat[[1=TS; 2=Correlation Mat's]][,,person]
X <- rbind( t(simdat[[1]][,,1]), t(simdat[[1]][,,2]), t(simdat[[1]][,,3]), t(simdat[[1]][,,4]) )
X
X <- array(0,  dim = c(N*tl, nts ) )

groupvec <- rep(c(1:N),  each = tl )

return_standat <- dcnet:::standat( data = X, J = N, group = groupvec, xC = groupvec,  P = 1,
                                  standardize_data = TRUE,
                                  Q = 1, distribution = 0, meanstructure = "VAR")

return_standat

stan_data <- return_standat[ c("T", "xC", "rts", "nt",
                                   "distribution", "P", "Q",
                                   "meanstructure", "J", "nobs", "group")]

stan_data
stan_data$rts
stan_data$rts <- array( apply(simdat[[1]], 1, FUN = rbind ), dim = c(tl, N, nts) )
## same as this from X:  array(X,  dim = c(20, 4, 3) )



simdat[[1]][,,4]

simdat[[1]]

#saveRDS(object =stan_data ,file = '../../TEST/stan_data')

parameterization <- "DCC"

stanmodel <- switch(parameterization,
                        CCC = dcnet:::stanmodels$CCCMGARCH,
                        DCC = dcnet:::stanmodels$DCCMGARCH,
                        NULL)

stanmodel

model_fit <- rstan::sampling(stanmodel,
                                     data = stan_data,
                                     verbose = TRUE,
                                     iter = 500,
                                     control = list(adapt_delta = .99),
                                     chains = 4,
                                     init_r = .05)


options(width = 220 )
#model_fit

print(model_fit, pars =  c( 'phi0_fixed')) ## What's up wit the stdnorms? Seem WAY to large



fit <- dcnet( data = X, J = 4, group = groupvec, iterations = 500, standardize_data = TRUE)

fit
