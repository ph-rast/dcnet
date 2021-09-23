library(dcnet )

## Create data:
simdat <- dcnet:::.simuDCC(tslength = 20,  N = 4,  n_ts = 3,  ranef_sd_S = 0.0001 )
simdat

## Generated TS for person 1: simdat[[1=TS; 2=Correlation Mat's]][,,person]
X <- rbind( t(simdat[[1]][,,1]), t(simdat[[1]][,,2]), t(simdat[[1]][,,3]), t(simdat[[1]][,,4]) )
X

X[,2] <- X[,2] + 10

groupvec <- rep(c(1:4),  each = 20 )

return_standat <- dcnet:::standat( data = X, J = 4, group = groupvec, xC = groupvec,  P = 1,
                                  standardize_data = TRUE,
                                  Q = 1, distribution = 0, meanstructure = "VAR")

return_standat

stan_data <- return_standat[ c("T", "xC", "rts", "nt",
                                   "distribution", "P", "Q",
                                   "meanstructure", "J", "nobs", "group")]

stan_data
stan_data$rts
stan_data$rts <-
    array(X,  dim = c(20, 4, 3) )

simdat[[1]][,,4]

simdat[[1]]

stan_data

parameterization <- "DCC"

stanmodel <- switch(parameterization,
                        CCC = dcnet:::stanmodels$CCCMGARCH,
                        DCC = dcnet:::stanmodels$DCCMGARCH,
                        NULL)

stanmodel

model_fit <- rstan::sampling(stanmodel,
                                     data = stan_data,
                                     verbose = TRUE,
                                     iter = 200,
                                     control = list(adapt_delta = .99),
                                     chains = 4,
                                     init_r = .05)

options(width = 180 )
model_fit
print(model_fit, pars =  'phi0')



fit <- dcnet( data = X, J = 4, group = groupvec, iterations = 500, standardize_data = TRUE)

fit
