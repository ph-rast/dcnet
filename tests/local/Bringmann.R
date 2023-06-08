## Data from Bringmann 2013: https://dx.plos.org/10.1371/journal.pone.0060188
## Reanalized by Epskamp https://arxiv.org/pdf/1609.04156.pdf
library(data.table )

# bman <- read.csv(file = "https://doi.org/10.1371/journal.pone.0060188.s004" )
# saveRDS(object = bman, file = "bman.Rds")
bman <- data.table( readRDS("bman.Rds") )
head(readRDS("bman.Rds"), n =  5)

bman

## informat04: From Bringmann's paper: "...equals to 0 if person belongs to the control group, and takes value 1 if the person received mindfulness therapy"
## keep complete cases and only those in  the treatment group
bmcomp <- bman[ complete.cases(bman )  &  bman$informat04 == 1]

length(unique(bmcomp$subjno))

## Select only individuals with 80+ entries
tmp <- bmcomp[, .N, by = subjno]
tmp
sel <- tmp[N >= 80]$subjno
length(sel)

dt <- bmcomp[subjno %in% sel]

## align so that everyone has same st_period transition
## Goal: Align everyone at the same 0/1 transition
## find shortest 0 sequence and truncate all others to that
short <- min( dt[, sum(st_period == 0), by = subjno]$V1 )
short
## Create starting var
dt[, start := sum(st_period == 0)-short+1, by = subjno]

## Create observation index
dt[, obs := 1:.N, by = subjno]
dt

dt[, diff := sign(obs-start)]
## keep only non-negatives
dt <- dt[dt$diff >= 0,]

## Sanit check: All zero length must be the same
dt[, sum(st_period == 0), by = subjno]$V1

tsl <- min( dt[, .N, by = subjno]$N )
tsl
## Truncate the data table to exactly tsl observations per person
## The .SD refers to the current group of rows being processed.
dt <- dt[, head(.SD, tsl), by = subjno]

## Re-Create observation index
dt[, obs := 1:.N, by = subjno]


library(ggplot2 )
ggplot(dt,  aes(x = obs,  y = opgewkt_, group = factor(subjno), color = factor(subjno)) ) + geom_smooth( show.legend = FALSE, se = FALSE) + geom_point(alpha = .2) + geom_vline(xintercept =  36 )
dev.off( )

head(dt )

N <- length(unique(dt$subjno ) )
N
ids <- unique(dt$subjno )
tl <- tsl
groupvec <- rep(seq_len(N),  each = tl)

varnames <- c("opgewkt_", "somber__", "ontspann" ,"onplplez")

tsdat <- lapply(ids, function(x) dt[dt$subjno == x, .SD, .SDcols = varnames] )
tsdat

scale(tsdat )

getwd( )
setwd( "../")

devtools::load_all( )

st_period <- lapply(ids, function(x) dt[dt$subjno == x, st_period] )
st_period[[1]][37]


## Model with just constnat cor
fit0 <- dcnet( data = tsdat, J = N, group = groupvec, S_pred = st_period, parameterization = "CCC",
             standardize_data = FALSE, sampling_algorithm = 'variational', threads = 1, init = 1)#, iterations = 10, threads = 1, init = 0)

## Intervention not modeled
fitN <- dcnet( data = tsdat, J = N, group = groupvec, S_pred = NULL, parameterization = "DCCr",
             standardize_data = FALSE, sampling_algorithm = 'variational', threads = 1, init = 1)#, iterations = 10, threads = 1, init = 0)

## Model with two S matrices pre/post intervention
fit <- dcnet( data = tsdat, J = N, group = groupvec, S_pred = st_period, parameterization = "DCCr",
             standardize_data = FALSE, sampling_algorithm = 'variational', threads = 1, init = 1)#, iterations = 10, threads = 1, init = 0)

nts <- length(varnames)
zeroord <- matrix(data.table(fit$model_fit$summary( variables = c("Sfixed")))$mean, ncol = nts)
-cov2cor( solve(zeroord ) )
data.table(fit$model_fit$summary( variables = c("phi0_fixed")))$mean

## Model Selection
## Extract log_lik
log_lik_loc0 <- grep( "log_lik", colnames(fit0$model_fit$draws( )) )
log_lik0 <- fit0$model_fit$draws( )[, log_lik_loc0]

log_lik_locN <- grep( "log_lik", colnames(fitN$model_fit$draws( )) )
log_likN <- fitN$model_fit$draws( )[, log_lik_locN]

log_lik_loc <- grep( "log_lik", colnames(fit$model_fit$draws( )) )
log_lik <- fit$model_fit$draws( )[, log_lik_loc]

f0 <- loo::loo( log_lik0)
fN <- loo::loo( log_likN)
fr <- loo::loo( log_lik)

loo::loo_compare(f0, fN,  fr)

## Plots
## R[id, timepoint, variable, variable]
## Take median of posterior distribution

ind <- seq_len(N)

## for HMC:
median(fit$model_fit$draws( )[,,'R[1,1,1,2]'])
## for VB
median(fit$model_fit$draws( )[,'R[1,1,1,2]'])

out <- NULL
for(i in ind) {
#  out <- c(out,sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[,, paste0('R[',i,',',x,',1,2]')])))
  out <- c(out,sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[,paste0('R[',i,',',x,',1,2]')])))
}
out

nts <- 4
## Nested loop through first row of correlations (1,2; 1,3; 1,4 ...) for all tl and all N
plout <- sapply(2:nts,  function(p) {
  sapply(ind,  function(f) {
    sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('R[',f,',',x,',1,',p,']')]))
  })
})

plout2 <- sapply(3:nts,  function(p) {
  sapply(ind,  function(f) {
    sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('R[',f,',',x,',2,',p,']')]))
  })
})

plout3 <- sapply(4,  function(p) {
  sapply(ind,  function(f) {
    sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('R[',f,',',x,',3,',p,']')]))
  })
})


plout <- sapply(2:nts,  function(p) {
  sapply(ind,  function(f) {
    sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('pcor[',f,',',x,',1,',p,']')]))
  })
})

plout2 <- sapply(3:nts,  function(p) {
  sapply(ind,  function(f) {
    sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('pcor[',f,',',x,',2,',p,']')]))
  })
})

plout3 <- sapply(4,  function(p) {
  sapply(ind,  function(f) {
    sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('pcor[',f,',',x,',3,',p,']')]))
  })
})


df <- data.frame( plout, time = rep(seq(1:tl), N ), id = as.factor(rep(ind,  each = tl)))
names(df)[1:3] <- c('cor12',  'cor13',  'cor14' )

names(df)[1:2] <- c('cor12',  'cor13' )
names(df)[1] <- c('cor12')

head(df)


library(ggplot2 )
varnames
#'enthus', 'fear', 'angry', 'happy'
c12 <-
  ggplot(df,  aes(x = time,  y = cor12 , color = id)) + geom_line(show.legend = FALSE )+ coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous( paste0("PCor(", varnames[1], ", ", varnames[2], ")") ) + geom_vline( xintercept = 36)
c13 <-
  ggplot(df,  aes(x = time,  y = cor13 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous( paste0("PCor(", varnames[1], ", ", varnames[3], ")") )+ geom_vline( xintercept = 36)
c14 <- ggplot(df,  aes(x = time,  y = cor14 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous( paste0("PCor(", varnames[1], ", ", varnames[4], ")") )+ geom_vline( xintercept = 36)



df2 <- data.frame( plout2, time = rep(seq(1:tl), N ), id = as.factor(rep(ind,  each = tl)))
names(df2)[1:2] <- c('cor23',  'cor24' )
head(df2)

c23 <- ggplot(df2,  aes(x = time,  y = cor23 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous( paste0("PCor(", varnames[2], ", ", varnames[3], ")") )+ geom_vline( xintercept = 36)
c24 <- ggplot(df2,  aes(x = time,  y = cor24 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous( paste0("PCor(", varnames[2], ", ", varnames[4], ")") )+ geom_vline( xintercept = 36)



df3 <- data.frame( plout3, time = rep(seq(1:tl), N ), id = as.factor(rep(ind,  each = tl)))
names(df3)[1] <- c('cor34')
head(df3)

c34 <- ggplot(df3,  aes(x = time,  y = cor34 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous(paste0("PCor(", varnames[3], ", ", varnames[4], ")"))+ geom_vline( xintercept = 36)

nn <- ggplot( ) + theme_void()

c12

ggplot(df3, aes(x = time,  y = cor34 , color = id)) + geom_line( )

library(patchwork )

(c12 | c13 | c14 ) /
( nn | c23 | c24 ) /
( nn | nn  | c34 )




### Yrep => rts_out[N,ts_length,timesries]
tsl
N
fit$model_fit$draws( )[,'rts_out[46,68,4]']

person <- 6
## extract lower, median, upper quantile and scal back to original scale
yrep4 <- sapply(seq_len(tl), function(x) {
  quantile(fit$model_fit$draws( )[, paste0('rts_out[',person,',',x,',1]')], c(.025, .5, .975))
  })


yrep4 <- data.frame(yrep = t(yrep4))
yrep4$time <- seq_len(tl )
names(yrep4)
yrep4$obs <- unlist( tsdat[[person]][,1] )

tsdat[[1]][,1]
str(yrep4 )

ggplot(yrep4,  aes(x = time, y = yrep.50.) ) + geom_line( ) + geom_ribbon(aes(ymin =  yrep.2.5.,  ymax = yrep.97.5.),  alpha = .2) + geom_line( aes(y = obs ), color = 'red')

## Sanity check to see if stan_model transform to original data
yrep4$standat <- fit$RTS_full[[person]][,1]
ggplot(yrep4,  aes(x = time, y = standat) ) + geom_point( ) + geom_ribbon(aes(ymin =  yrep.2.5.,  ymax = yrep.97.5.),  alpha = .2) + geom_line( aes(y = obs ), color = 'red')




fit$RTS_full


### Check ppd

y_rep_loc <- grep( "rts_out", colnames(fit$model_fit$draws( )) )
y_rep <-fit$model_fit$draws( )[, y_rep_loc]
str(y_rep)

plot(density(as.numeric( y_rep[1, ])), col = '#33333375', ylim = c(0, 0.6))
for(i in sample(nrow(y_rep), size = 20)){
    lines(density( as.numeric(y_rep[i,]) ) , col = '#33333375', lwd = 1)
}

lines(density(unlist(fit$RTS_full)), col = 'red', lwd = 3)

dev.off()



####
library(qgraph )

## Plot unconditional Corrs:
zeroord <- matrix(data.table(fit$model_fit$summary( variables = c("Sfixed")))$mean, ncol = nts)
ucpcor <- -cov2cor( solve(zeroord ) )
diag(ucpcor ) <- -diag(ucpcor )

qgraph::qgraph(zeroord, labels = varnames)
qgraph::qgraph(ucpcor, labels = varnames)


## individual, time-varying:
head(df )
head(df2)
head(df3)


rmats <- function(i ){
  indR <- diag(1, nrow = 4 )
  indR[upper.tri(indR)] <-  c( unlist(df[i,1:3]), unlist(df2[i, 1:2]), unlist(df3[i, 1]) )
  indR <- Matrix::forceSymmetric(indR )
  return(indR )
}

tsl
indRlist <- lapply(1:tsl, FUN = rmats)

## needs to be list of full cormats
qgraph::qgraph.animate( input = indRlist )


for(k in seq_len(tsl) ){
  png(filename = paste0("~/Downloads/anim/abc_", sprintf("%02d", k), ".png"))
  qgraph::qgraph(indRlist[[k]], labels = varnames)
  dev.off( ) 
}

#system( command = "convert -delay 1 -loop 0 ~/Downloads/anim/abc* anim.gif" )


### VAR part

grep('phi0', colnames(fit$model_fit$draws( )) )



colMeans( fit$model_fit$draws( )[, 4552:4600] )
colMeans( fit$model_fit$draws( )[, 1057:1072] )

mS <- matrix( colMeans(  fit$model_fit$draws( )[, 1041:1056] ), ncol =  4, byrow =  TRUE)
mS2 <- matrix( colMeans( fit$model_fit$draws( )[, 1057:1072] ), ncol =  4, byrow =  TRUE)

mS2

fit$model_fit$draws( )[, c( 'S[1,2]','S[1,3]','S[1,4]','S[2,3]','S[2,4]','S[3,4]')]


qgraph::qgraph(mS2, labels = varnames)






####
x <- tanh(rnorm(6,0,1) + rgamma(6,  1,  1 )*rnorm(6, 0, 1) )
x <- tanh(rnorm(6,0,1) + rgamma(6, 0.05, 2) )
nt <- 4
choose(4, 2)
z <- diag(rep(0, nt))
z

index <- 0
for(j in 1:nt) {
    for(i in 1:nt) {
      if( i > j ) {
	index = index + 1
	z[i,j] = x[index] 
      }
      if (i < j) {
	z[i,j] = 0
	  }
    }
}
z

L <- diag(rep(0, nt))
for(i in 1:nt) {
    for(j in 1:nt){
      if(i < j){
	L[i,j]=0.0;
      }
      if(i == j){
	if(i == 1){
	  L[i,j]=1.0;
	    }
	if(i > 1){ 
	  L[i,j]=sqrt( 1 - sum( L[i, 1:j]^2 ) );
	}
      }
      if(i > j){
	L[i,j]=z[i,j]*sqrt(1 - sum( L[i, 1:j]^2) );
      }
    }
}
L
R <- L%*%t(L)
round(R, 3)
