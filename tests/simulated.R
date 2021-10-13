library(dcnet )

## Create data:
N <- 20
tl <- 40
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
                                     iter = 100,
                                     #control = list(adapt_delta = .99),
                                     chains = 1,
                                     init_r = .05)


options(width = 220 )
#model_fit

print(model_fit, pars =  c("H")) ## What's up wit the stdnorms? Seem WAY to large

tail( rstan::summary(model_fit, pars = "H")$summary , n =  18)
N
H <- array( rstan::summary(model_fit, pars = "H")$summary[,1], dim = c(N, tl, 3, 3 ) )
H <- array( rstan::summary(model_fit, pars = "H")$summary[,"mean"], dim = c(3, 3, N, tl) )

r1 <- array(NA, dim = c(N, tl ) )
r2 <- r1
r3 <- r1

for( i in 1:N ) {
    for(k in 1:tl ) {
        r1[i,k] <- cov2cor(H[,,i,k])[2,1]
        r2[i,k] <- cov2cor(H[,,i,k])[3,1]
        r3[i,k] <- cov2cor(H[,,i,k])[3,2]
    }
}

plot(1:tl, r1[1,], type =  'l', ylim = c(-1, 1 ), col = "#3299ff60")
for(i in 2:N ) {
    lines(r1[i,],  col = "#3299ff60" )
}


op <- par(mfcol = c(2, 2 ) )
plot(1:tl, r1[1,], type =  'l', ylim = c(-1, 1 ), col = 2)
for(i in 2:N ) {
    lines(r1[i,],  col = i)#paste0("#99f", sprintf("%02d", 1), "f60") )
}

plot(1:tl, r2[1,], type =  'l', ylim = c(-1, 1 ), col = 2)
for(i in 2:N ) {
    lines(r2[i,],  col = i)#paste0("#99f", sprintf("%02d", 1), "f60") )
}

plot(1:tl, r3[1,], type =  'l', ylim = c(-1, 1 ), col = 2)
for(i in 2:N ) {
    lines(r3[i,],  col = i)#paste0("#99f", sprintf("%02d", 1), "f60") )
}

op


mnr1 <-NA# array(NA, dim = c(1,  tl ) )
## Average over N's
for( k in 1:tl ) {
    mnr1[k] <- mean(H[1,2,,k]) /( sqrt( mean(H[1,1,,k]) * mean(H[2,2,,k]) ) )
}
mnr1

plot(1:tl, r1[1,], type =  'l', ylim = c(-1, 1 ), col = "#5033ff60")
for(i in 2:N ) {
    lines(r1[i,],  col = "#5033ff60") 
}
lines(mnr1, col = 1)

print(model_fit, pars =  c("S")) 



fit <- dcnet( data = X, J = 4, group = groupvec, iterations = 500, standardize_data = TRUE)

fit
