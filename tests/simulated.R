#library(dcnet )

## Create data:
N <- 5
tl <- 50
nts <- 3
simdat <- dcnet:::.simuDCC(tslength = tl,  N = N,  n_ts = nts,  ranef_sd_S = 0.0001, phi0_fixed =  c(1, 10, 20 ))
dim(simdat[[1]])

rtsgen <- lapply(seq(dim(simdat[[1]])[3]), function(x) t( simdat[[1]][,,x] ))

## Generated TS for person 1: simdat[[1=TS; 2=Correlation Mat's]][,,person]
X <- rbind( t(simdat[[1]][,,1]), t(simdat[[1]][,,2]), t(simdat[[1]][,,3]), t(simdat[[1]][,,4]) )
X
X <- array(0,  dim = c(N*tl, nts ) )

groupvec <- rep(c(1:N),  each = tl )
groupvec

return_standat <- dcnet:::standat( data = X, J = N, group = groupvec, xC = groupvec,  P = 1,
                                  standardize_data = TRUE,
                                  Q = 1, distribution = 0, meanstructure = "VAR")

return_standat

stan_data <- return_standat[ c("T", "xC", "rts", "nt",
                                   "distribution", "P", "Q",
                                   "meanstructure", "J", "nobs", "group")]

stan_data
stan_data$rts
stan_data$rts <- rtsgen

## same as this from X:  array(X,  dim = c(20, 4, 3) )



stan_data$rts


#saveRDS(object =stan_data ,file = '../../TEST/stan_data')

parameterization <- "VAR"

stanmodel <- switch(parameterization,
                    DCCf = dcnet:::stanmodels$DCCMGARCHfixedS,
                    DCC = dcnet:::stanmodels$DCCMGARCH,
                    VAR = dcnet:::stanmodels$VAR,
                    NULL)

stanmodel


system.time( {
  model_fit <- rstan::sampling(stanmodel,
                             data = stan_data,
                             verbose = TRUE,
                             iter = 200,
                                        #warmup =  1500, 
                             control = list(adapt_delta = .99),
                             chains = 4,
                             init_r = .05)
})

library(cmdstanr )
file <- file.path("../inst/stan/VAR.stan" )
mod <- cmdstan_model(file)

model_fit <- mod$sample(
                   data = stan_data,
                   chains = 4,
                   parallel_chains = 4,
                   iter_warmup = 500,
                   iter_sampling = 500)

                                        #model_fit <- rstan::vb(stanmodel, data = stan_data,
                       ##algorithm = "fullrank",
                       #iter = 20000, init_r = 0.01)


options(width = 220 )
model_fit$summary( )
model_fit$summary("phi" )


print(model_fit, pars =  c("rescov"))
print(model_fit, pars =  c("phi"))
print(model_fit, pars =  c("phi0_fixed"))
print(model_fit, pars =  c("phi0_tau"))
print(model_fit, pars =  c("vec_phi_fixed"))

print(model_fit, pars =  c("H")) ## What's up wit the stdnorms? Seem WAY to large
print(model_fit, pars =  c("S_Lv_fixed")) ## What's up wit the stdnorms? Seem WAY to large
print(model_fit, pars =  c("S")) ## What's up wit the stdnorms? Seem WAY to large


slv <- rstan::summary(model_fit, pars = "S_Lv_fixed")$summary[,'97.5%']
slv
slv <- c(1,  slv )
slv
L <- matrix(0, ncol = nts, nrow = nts)
L

index <- 0
for( j in 1:nts ) {
    for( i in 1:nts ) {
        if(i >=  j) {
            index <-  index + 1
            L[i,j] <- slv[index]
        } else next
    }
}

cov2cor( L %*% t(L) )

                                        #H <- array( rstan::summary(model_fit, pars = "H")$summary[,1], dim = c(N, tl, 3, 3 ) )
nts
H <- array( rstan::summary(model_fit, pars = "H")$summary[,"mean"], dim = c(nts, nts, N, tl) )

r21 <- array(NA, dim = c(N, tl ) )
r31 <- r21
r41 <- r21
r32 <- r21
r42 <- r21
r43 <- r21


for( i in 1:N ) {
    for(k in 1:tl ) {
        r21[i,k] <- cov2cor(H[,,i,k])[2,1]
        r31[i,k] <- cov2cor(H[,,i,k])[3,1]
        r41[i,k] <- cov2cor(H[,,i,k])[4,1]
        r32[i,k] <- cov2cor(H[,,i,k])[3,2]
        r42[i,k] <- cov2cor(H[,,i,k])[4,2]
        r43[i,k] <- cov2cor(H[,,i,k])[4,3]
    }
}

r21
#Sys.setenv("DISPLAY"=":0.0")

pdf(file = "./local/corr.pdf", width = 10, height =  6)

plot(1:tl, r21[1,], type =  'l', ylim = c(-1, 1 ), col = "#3299ff60")
for(i in 2:N ) {
    lines(r21[i,],  col = "#3299ff60" )
}

dev.off( )


pdf(file = "./local/corr.pdf", width = 10, height =  6)

op <- par(mfcol = c(2, 3 ) )
plot(1:tl, r21[1,], type =  'l', ylim = c(-1, 1 ), col = 2)
for(i in 2:N ) {
    lines(r21[i,],  col = i)#paste0("#99f", sprintf("%02d", 1), "f60") )
}
plot(1:tl, r31[1,], type =  'l', ylim = c(-1, 1 ), col = 2)
for(i in 2:N ) {
    lines(r31[i,],  col = i)#paste0("#99f", sprintf("%02d", 1), "f60") )
}
plot(1:tl, r41[1,], type =  'l', ylim = c(-1, 1 ), col = 2)
for(i in 2:N ) {
    lines(r41[i,],  col = i)#paste0("#99f", sprintf("%02d", 1), "f60") )
}
plot(1:tl, r32[1,], type =  'l', ylim = c(-1, 1 ), col = 2)
for(i in 2:N ) {
    lines(r32[i,],  col = i)#paste0("#99f", sprintf("%02d", 1), "f60") )
}
plot(1:tl, r42[1,], type =  'l', ylim = c(-1, 1 ), col = 2)
for(i in 2:N ) {
    lines(r42[i,],  col = i)#paste0("#99f", sprintf("%02d", 1), "f60") )
}
plot(1:tl, r43[1,], type =  'l', ylim = c(-1, 1 ), col = 2)
for(i in 2:N ) {
    lines(r43[i,],  col = i)#paste0("#99f", sprintf("%02d", 1), "f60") )
}
op

dev.off( )



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





## Real data: Nestler
dfs <- read.table("./local/dataFlip.txt", header = TRUE ) 
head( dfs )

dat1 <- dfs[,c("id","day", "gm_sociable",  "gm_creative", "gm_friendly", "gm_organised")]

## make wide to create NA's and equal length for all
X0 <- reshape(dat1, idvar = "id", timevar = "day", direction = "wide",
        v.names = c("gm_sociable",  "gm_creative", "gm_friendly", "gm_organised"))

## REplace missings with 0
X0[is.na(X0 )] <- 0

head(X0 )
length(names(X0))

XL <- reshape(X0, direction = "long",
              idvar = "id", sep = ".",
              varying = names(X0)[2:329] )

## sort by id, then time
X2 <- XL[order( XL$id,  XL$time), ]
head(X2 )
unique(X2$id )
X2 <- X2[X2$id <= 50,]
X2 <- X2[X2$time <= 24,]

N <- max( X2$id )
N
groupvec <- X2$id
groupvec
tl <- dim(X2 )[1]/ N
tl

length(groupvec )
nts <- 4 # Number of variables

return_standat <- dcnet:::standat( data = X2[,c("gm_sociable",  "gm_creative", "gm_friendly", "gm_organised")],
                                  J = N, group = groupvec, xC = groupvec,  P = 1,
                                  standardize_data = TRUE,
                                  Q = 1, distribution = 0, meanstructure = "VAR")


stan_data <- return_standat[ c("T", "xC", "rts", "nt",
                                   "distribution", "P", "Q",
                                   "meanstructure", "J", "nobs", "group")]



rts <- array(NA, dim = c(tl, N, nts))


rts[,,1] <- array( X2[,c("gm_sociable")], dim = c(tl, N) ) ## time by person
rts[,,2] <- array( X2[,c("gm_creative")], dim = c(tl, N) ) ## time by person
rts[,,3] <- array( X2[,c("gm_friendly")], dim = c(tl, N) ) ## time by person
rts[,,4] <- array( X2[,c("gm_organised")], dim = c(tl, N) ) ## time by person

stan_data$rts <- rts

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
                                     control = list(adapt_delta = .95),
                                     chains = 4)#,
                                     #init_r = .05)


#model_fit <- rstan::vb(stanmodel, data = stan_data,
#                             iter =  10000)


options(width = 220 )
#model_fit
