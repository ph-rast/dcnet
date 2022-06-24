##library(dcnet )
devtools::load_all( )
options(width = 200 )

## Create data:
N <- 50
tl <- 100
nts <- 4
simdat <- .simuDCC(tslength = tl,  N = N,  n_ts = nts,  ranef_sd_S = 0.1,
                   phi0_fixed =  c(0, 0, 0 , 0), alpha = .5)

rtsgen <- lapply(seq(dim(simdat[[1]])[3]), function(x) t( simdat[[1]][,,x] ))
## Fixed Corr
simdat$S

typeof(rtsgen )
dim(rtsgen[[1]])
rtsgen[[1]]

## Generated TS for person 1: simdat[[1=TS; 2=Correlation Mat's]][,,person]
X <- rbind( t(simdat[[1]][,,1]), t(simdat[[1]][,,2]), t(simdat[[1]][,,3]), t(simdat[[1]][,,4]) )

X <- array(0,  dim = c(N*tl, nts ) )

groupvec <- rep(c(1:N),  each = tl )


fit <- dcnet( data =  rtsgen, J =  N, group =  groupvec, standardize_data = FALSE,
             iterations = 400, sampling_algorithm = 'HMC')

dcnet::.get_stan_summary(fit$model_fit, )

devtools::load_all( )
summary(fit )


dcnet::print.summary.dcnet(fit )
print.summary.dcnet(out)






###################
## Manual testing

fit$model_fit$draws( )
fit$sampling_algorithm


summary.dcnet(fit)

fit$sampling_algorithm


print.summary.dcnet( )


fit$TS_names


 metaNames <- c("param", "distribution", "num_dist", "iter", "chains",
                   "elapsed_time", "date", "nt", "TS_names", "mgarchQ",
                   "mgarchP", "meanstructure", "sampling_algorithm")

with(fit,  mget(metaNames) )

fit$param
object$sampling_algorithm
sum(grepl("CmdStanFit", class(fit$model_fit )))


rtsgen
stan_data <- stan_data( data =  rtsgen, J =  N, group =  groupvec, xC = NULL, standardize_data = FALSE)
colnames(stan_data$rts)

## fit is a tibble object
tbl <- fit$model_fit$summary( mn =  ~ mean(.), sd =  ~ sd(.), qtl = ~ quantile(., probs = CrI ))
tbl

fit$model_fit$draws( )
fit$model_fit$time( )


#library(tidyverse )
#tbl %>% filter(grepl('lp|phi0', variable))


## Go with data.table
library(data.table )
dt <- as.data.table(tbl)

## cf. https://www.geeksforgeeks.org/select-rows-with-partial-string-match-in-r-dataframe/
## common and arma params
ap <- "lp|phi0_fixed|phi0_L|phi0_tau|vec_phi_fixed|phi_L|phi_tau"

## DCC params:
## D
dccD <- "c_h_fixed|a_h_fixed|b_h_fixed|c_h_L|c_h_tau|a_h_L|a_h_tau|b_h_L|b_h_tau"
## Q
dccQ <- "l_a_q|l_a_q_r|l_b_q|l_b_q_r|S"

paste0(dccD, '|', dccQ)

dt[ variable %like% dccQ ]







library(cmdstanr )
getwd( )
file <- file.path("../inst/stan/DCCMGARCHfixedS.stan" ) ## convergence with variational
file <- file.path("../inst/stan/DCCMGARCHrandQ.stan" )
##file <- file.path("../inst/stan/DCCMGARCHfixedD.stan" ) ##
file

mod <- cmdstan_model(file, include_paths = "../inst/stan/")

model_fit <-mod$optimize( data = return_standat)

model_fit <-mod$variational(
                  data = return_standat,
##                   threads = 6,
##                  tol_rel_obj = 0.005,
                   iter =  30000)

model_fit$output( )

model_fit <- mod$sample(
                   data = stan_data,
                   chains = 4,
                   parallel_chains = 4,
                   iter_warmup = 1000,
                   iter_sampling = 1000,
                   adapt_delta = 0.95)


options(width = 220 )

model_fit$summary("l_a_q" )
model_fit$summary("l_b_q" )
model_fit$summary("l_a_q_r" )
model_fit$summary("a_q" )
model_fit$summary("b_q" )

model_fit$summary("phi" )
model_fit$print("phi0_fixed", max_rows = 120)

model_fit$print("S", max_rows = 90)

model_fit$print("H", max_rows = 200)
model_fit$print("R", max_rows = 500)

## individual draws.. not necessary I think
model_fit$draws()[,, 1]
lapply(1:10,  function(x) model_fit$draws()[,, x] )

which(unlist(model_fit$summary()[1]) == "R[1,1,1,1]")
df <- as.data.frame(model_fit$summary())
head(df)

which(df$variable == "R[1,1,2,1]")
df[2655:2700,]

model_fit$summary("c_h_fixed" )
model_fit$summary("c_h_tau" )
model_fit$summary("c_h_stdnorm" )



print(model_fit, pars =  c("rescov"))
print(model_fit, pars =  c("phi"))
print(model_fit, pars =  c("phi0_fixed"))
print(model_fit, pars =  c("phi0_tau"))
print(model_fit, pars =  c("vec_phi_fixed"))

print(model_fit, pars =  c("H")) ## What's up wit the stdnorms? Seem WAY to large
print(model_fit, pars =  c("S_Lv_fixed")) ## What's up wit the stdnorms? Seem WAY to large
print(model_fit, pars =  c("S")) ## What's up wit the stdnorms? Seem WAY to large


print(model_fit, pars =  c("H")) ## What's up wit the stdnorms? Seem WAY to large

draws <- model_fit$draws("R")
df <- data.frame(posterior::as_draws_df(draws))


## Person, Time, var, var
quantile(df$R.1.1.2.1., c(.05, .5, .95))
quantile(df$R.1.2.2.1., c(.05, .5, .95))
quantile(df$R.1.3.2.1., c(.05, .5, .95))
quantile(df$R.1.4.2.1., c(.05, .5, .95))



quantile(df$R.1.1.2.1., c(.05, .5, .95))




#slv <- rstan::summary(model_fit, pars = "S_Lv_fixed")$summary[,'97.5%']
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
#H <- array( rstan::summary(model_fit, pars = "H")$summary[,"mean"], dim = c(nts, nts, N, tl) )

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
devtools::load_all( )

dfs <- read.table("./local/dataFlip.txt", header = TRUE ) 
head( dfs )

dat1 <- dfs[,c("id","day", "gm_sociable",  "gm_creative", "gm_friendly", "gm_organised")]

## make wide to create NA's and equal length for all
X0 <- reshape(dat1, idvar = "id", timevar = "day", direction = "wide",
        v.names = c("gm_sociable",  "gm_creative", "gm_friendly", "gm_organised"))

X0
## Replace missings with 0
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
X2$time
X2 <- X2[X2$id <= 50,]
X2 <- X2[X2$time <= 70,]

N <- max( X2$id )
N
groupvec <- X2$id
groupvec
tl <- dim(X2 )[1]/ N
tl

length(groupvec )
nts <- 4 # Number of variables

## X2 needs to be a list of N matrices with dimension ntsXtl 
## Drop id and time variable
tsdat <- lapply( seq_len(N), function(x) X2[X2$id == x, 3:6])
str(tsdat)


fit <- dcnet( data = tsdat, J =  N, group =  groupvec, standardize_data = TRUE,
             sampling_algorithm = 'variational')

summary(fit)

getwd( )
saveRDS( fit,  file = './local/fit_nestler.Rds')

## Plots
## R[id, timepoint, variable, variable]
## Take median of posterior distribution

out <- NULL
for(i in 1:N) {
  out <- c(out,sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[,paste0('R[',i,',',x,',1,2]')])))
}


## Nested loop through first row of correlations (1,2; 1,3; 1,4 ...) for all tl and all N
plout <- sapply(2:nts,  function(p) {
  sapply(1:N,  function(f) {
    sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('R[',f,',',x,',1,',p,']')]))
  })
})


df <- data.frame( plout, time = rep(seq(1:tl), N ), id = as.factor(rep(1:N,  each = tl)))
names(df)[1:3] <- c('cor12',  'cor13',  'cor14' )
head(df)


library(ggplot2 )

c12 <- ggplot(df,  aes(x = time,  y = cor12 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +  scale_y_continuous("PCor(Sociable, Creative)")
c13 <- ggplot(df,  aes(x = time,  y = cor13 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +  scale_y_continuous("PCor(Sociable, Friendly)")
c14 <- ggplot(df,  aes(x = time,  y = cor14 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) )+  scale_y_continuous("PCor(Sociable, Organised)")


plout2 <- sapply(3:nts,  function(p) {
  sapply(1:N,  function(f) {
    sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('R[',f,',',x,',2,',p,']')]))
  })
})

df2 <- data.frame( plout2, time = rep(seq(1:tl), N ), id = as.factor(rep(1:N,  each = tl)))
names(df2)[1:2] <- c('cor23',  'cor24' )
head(df2)

c23 <- ggplot(df2,  aes(x = time,  y = cor23 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +  scale_y_continuous("PCor(Creative, Friendly)")
c24 <- ggplot(df2,  aes(x = time,  y = cor24 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +  scale_y_continuous("PCor(Creative, Organised)")


plout3 <- sapply(4,  function(p) {
  sapply(1:N,  function(f) {
    sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('R[',f,',',x,',3,',p,']')]))
  })
})

df3 <- data.frame( plout3, time = rep(seq(1:tl), N ), id = as.factor(rep(1:N,  each = tl)))
names(df3)[1] <- c('cor34')
head(df3)

c34 <- ggplot(df3,  aes(x = time,  y = cor34 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) )+  scale_y_continuous("PCor(Friendly, Organised)")

nn <- ggplot( ) + theme_void()

library(patchwork )

(c12 | c13 | c14 ) /
( nn | c23 | c24 ) /
( nn | nn  | c34 )  




## OLD
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


## model_fit <- rstan::sampling(stanmodel,
##                                      data = stan_data,
##                                      verbose = TRUE,
##                                      iter = 100,
##                                      control = list(adapt_delta = .95),
##                                      chains = 4)#,
##                                      #init_r = .05)


#model_fit <- rstan::vb(stanmodel, data = stan_data,
#                             iter =  10000)


options(width = 220 )
#model_fit
