## Data from: https://osf.io/emwr9/files/
## Howe, Reeves, & Fisher
options(width = 200, pillar.subtle = FALSE)

ptsd <- read.csv( file = "./2022_03_18_MCFA_data.csv" )
head( ptsd )

ptsd_comp <- ptsd#[complete.cases( ptsd ), ]
head( ptsd_comp )

## Select some variables
dat1 <- ptsd_comp[, c('id', 'enthus', 'fear', 'angry', 'happy',
                      'positive', 'horror', 'agg', 'shame',
                      'calm', 'morning','eve','night', 'sleepy', 'fatigue', 'mem')]
head(dat1 )
dat1[dat1$id == "P002",]

## drop all NA rowwise
dim(dat1 )
dat1c <- dat1[complete.cases(dat1),]

## add time variable to index all timepoints
ids <- unique( dat1c$id )
dat1c$time <- NA

for( i in 1:length(ids) ) {
  p_obs <- nrow( dat1c[dat1c$id == ids[i], ] )
  dat1c[dat1c$id == ids[i], 'time'] <- seq_len(p_obs)
}

## make wide to create NA's and equal length for all
X0 <- reshape(dat1c, idvar = "id", timevar = "time", direction = "wide",
              v.names = c('morning','eve', 'night','enthus', 'fear', 'angry', 'happy', 'positive', 'horror', 'agg', 'shame', 'calm','sleepy', 'fatigue', 'mem'))
head(X0)
## Replace missings with 50
#X0[is.na(X0 )] <- 50



X0$id <- seq_len(nrow(X0))

XL <- reshape(X0, direction = "long",
              idvar = "id", sep = ".",
              varying = names(X0)[2:length(names(X0))] )

XL
## sort by id, then time
X2 <- XL[order( XL$id,  XL$time), ]
head(X2, n =  10)

##
X2 <- X2[X2$time <= 95,]
sum(is.na(X2))
dim(X2)[1]

## Replace all missing with 50...  lol
#X2[is.na(X2 )] <- 50


N <- max( X2$id )
N

groupvec <- X2$id
groupvec

tl <- dim(X2 )[1]/ N
tl

nts <- 4 # Number of variables

## X2 needs to be a list of N matrices with dimension ntsXtl 
## Drop id and time variable
c('enthus', 'fear', 'angry', 'happy', 'positive', 'horror', 'agg', 'shame', 'calm')
varnames <- c('horror', 'sleepy', 'fatigue', 'mem')
names(X2 )
tsdat <- lapply( seq_len(N), function(x) X2[X2$id == x, varnames])
str(tsdat)

tsdat

devtools::load_all( )

X2$evenight <- ifelse(X2[, c('eve')] == 1 | X2[, c('night')] == 1, 1, 0)

S_pred <- lapply(seq_len(N), function(x) X2[X2$id == x, 'morning'] )
S_pred <- lapply(seq_len(N), function(x) c(rep(0,  47), rep(1, 48 )) )
S_pred <- lapply(seq_len(N), function(x) X2[X2$id == x, 'evenight'] )

S_pred

getwd( )
setwd( "../")
tsdat

fit <- dcnet( data = tsdat, J =  N, group =  groupvec, S_pred = NULL,
             standardize_data = TRUE, sampling_algorithm = 'variational')

summary(fit)

## Diff in S and S2

Ds12 <- fit$model_fit$draws( )[, 'S[1,2]'] - fit$model_fit$draws( )[, 'S2[1,2]']
Ds13 <- fit$model_fit$draws( )[, 'S[1,3]'] - fit$model_fit$draws( )[, 'S2[1,3]']
Ds14 <- fit$model_fit$draws( )[, 'S[1,4]'] - fit$model_fit$draws( )[, 'S2[1,4]']
Ds23 <- fit$model_fit$draws( )[, 'S[2,3]'] - fit$model_fit$draws( )[, 'S2[2,3]']
Ds24 <- fit$model_fit$draws( )[, 'S[2,4]'] - fit$model_fit$draws( )[, 'S2[2,4]']
Ds34 <- fit$model_fit$draws( )[, 'S[3,4]'] - fit$model_fit$draws( )[, 'S2[3,4]']

plot(density(Ds14 ) )
quantile( Ds12, probs = c(.05, .95) )
quantile( Ds13, probs = c(.05, .95) )
quantile( Ds14, probs = c(.05, .95) )
quantile( Ds23, probs = c(.05, .95) )
quantile( Ds24, probs = c(.05, .95) )
quantile( Ds34, probs = c(.05, .95) )




## Extract log_lik
log_lik_loc <- grep( "log_lik", colnames(fit$model_fit$draws( )) )
log_lik <- fit$model_fit$draws( )[, log_lik_loc]

## Gelman et al. BD3 p. 169 and pp. 173
waic <- function(log_lik ) {
  S <- nrow(log_lik)
  elpd <- sum( log( colSums(exp( log_lik )) / S ) )
  pwaic2 <- sum( apply( log_lik, 2, var) )  
  waic <- -2 * (elpd - pwaic2)
  return( list(waic = waic,  elpd = elpd, pwaic2 = pwaic2) )
}

waic(log_lik = log_lik)
loo::waic( log_lik)
loo::loo( log_lik)

fit_const <- dcnet( data = tsdat, J =  N, group =  groupvec, parameterization = 'CCC',
                   standardize_data = FALSE,
                   sampling_algorithm = 'variational')


log_lik_loc <- grep( "log_lik", colnames(fit_const$model_fit$draws( )) )
log_lik_const <- fit_const$model_fit$draws( )[, log_lik_loc]
waic( log_lik = fit_const$model_fit$draws( )[, log_lik_loc] )

loo::waic( log_lik_const )
m1 <- loo::loo( log_lik )
m2 <- loo::loo( x = log_lik_const)

loo::loo_compare(m1,  m2 )

rts_out_loc <- grep( "rts_out", colnames(fit$model_fit$draws( )) )

rts_out <- fit$model_fit$draws( )[, rts_out_loc]
head(rts_out)
## Format of variable names: rts_out[j, T, nt ]
rts_out[, 18:25]

#saveRDS( fit,  file = './local/ptsd.Rds')
#load(file = './local/ptsd.Rds')

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

plout2 <- sapply(3:nts,  function(p) {
  sapply(1:N,  function(f) {
    sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('R[',f,',',x,',2,',p,']')]))
  })
})

plout3 <- sapply(4,  function(p) {
  sapply(1:N,  function(f) {
    sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('R[',f,',',x,',3,',p,']')]))
  })
})


plout <- sapply(2:nts,  function(p) {
  sapply(1:N,  function(f) {
    sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('pcor[',f,',',x,',1,',p,']')]))
  })
})

plout2 <- sapply(3:nts,  function(p) {
  sapply(1:N,  function(f) {
    sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('pcor[',f,',',x,',2,',p,']')]))
  })
})

plout3 <- sapply(4,  function(p) {
  sapply(1:N,  function(f) {
    sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('pcor[',f,',',x,',3,',p,']')]))
  })
})




df <- data.frame( plout, time = rep(seq(1:tl), N ), id = as.factor(rep(1:N,  each = tl)))
names(df)[1:3] <- c('cor12',  'cor13',  'cor14' )
head(df)


library(ggplot2 )
varnames
#'enthus', 'fear', 'angry', 'happy'
c12 <- ggplot(df,  aes(x = time,  y = cor12 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous( paste0("PCor(", varnames[1], ", ", varnames[2], ")") )
c13 <- ggplot(df,  aes(x = time,  y = cor13 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous( paste0("PCor(", varnames[1], ", ", varnames[3], ")") )
c14 <- ggplot(df,  aes(x = time,  y = cor14 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous( paste0("PCor(", varnames[1], ", ", varnames[4], ")") )



df2 <- data.frame( plout2, time = rep(seq(1:tl), N ), id = as.factor(rep(1:N,  each = tl)))
names(df2)[1:2] <- c('cor23',  'cor24' )
head(df2)

c23 <- ggplot(df2,  aes(x = time,  y = cor23 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous( paste0("PCor(", varnames[2], ", ", varnames[3], ")") )
c24 <- ggplot(df2,  aes(x = time,  y = cor24 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous( paste0("PCor(", varnames[2], ", ", varnames[4], ")") )



df3 <- data.frame( plout3, time = rep(seq(1:tl), N ), id = as.factor(rep(1:N,  each = tl)))
names(df3)[1] <- c('cor34')
head(df3)

c34 <- ggplot(df3,  aes(x = time,  y = cor34 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous(paste0("PCor(", varnames[3], ", ", varnames[4], ")"))

nn <- ggplot( ) + theme_void()

ggsave(filename = "/home/philippe/UZH/Kongresse/Kongresse2022/Talk/Figures/pcor12.pdf",  plot = c12,  width = 5, height = 3 )


ggplot(df3, aes(x = time,  y = cor34 , color = id)) + geom_line( )

library(patchwork )

(c12 | c13 | c14 ) /
( nn | c23 | c24 ) /
  ( nn | nn  | c34 )

ggsave(filename = "/home/philippe/UZH/Kongresse/Kongresse2022/Talk/Figures/pcor.pdf", width = 8, height = 5.5)




## Compare Density


y_rep_loc <- grep( "rts_out", colnames(fit$model_fit$draws( )) )
y_rep <-fit$model_fit$draws( )[, y_rep_loc]
str(y_rep)

plot(density(as.numeric( y_rep[1, ])), col = '#33333375', ylim = c(0, 0.6))
for(i in sample(nrow(y_rep), size = 20)){
    lines(density( as.numeric(y_rep[i,]) ) , col = '#33333375', lwd = 1)
}

lines(density(unlist(fit$RTS_full)), col = 'red', lwd = 3)

dev.off()


## regexp variable 1
sel <- grep( pattern = "rts_out\\[[0-9]+,[0-95]+,1" , x = colnames(y_rep), perl = TRUE)
dim(y_rep[,sel])

y_rep_1 <- data.table( cbind(yrep = colMeans( y_rep[,sel] ), time = rep(1:95, N ), id = rep(1:N, each = 95 ) ))
y_rep_1$id <- as.factor(y_rep_1$id )
head(y_rep_1)

ggplot(y_rep_1,  aes(x = time, y = yrep, color = id)) + geom_line( )


id <- 4
df2 <- data.table(out = fit_const$RTS_full[[id]][, "horror"], time = seq(1:95))
head(df2 )
ggplot(y_rep_1[y_rep_1$id == paste0(id),],  aes(x = time, y = yrep)) + geom_line( ) + geom_line( data = df2, aes(y = out), color =  'red')
## CHECK hwo if RTS if manipulated in dcnet function

df2 <- data.table(out = tsdat[[1]][, "horror"], time = seq(1:95))
ggplot(y_rep_1[y_rep_1$id == "1",],  aes(x = time, y = yrep)) + geom_line( ) + geom_line( data = df2, aes(y = out), color =  'red')
## CHECK hwo if RTS if manipulated in dcnet function


library(qgraph )

grep('S', colnames(fit$model_fit$draws( )) )
colMeans( fit$model_fit$draws( )[, 1041:1056] )
colMeans( fit$model_fit$draws( )[, 1057:1072] )

mS <- matrix( colMeans( fit$model_fit$draws( )[, 1041:1056] ), ncol =  4, byrow =  TRUE)
mS2 <- matrix( colMeans( fit$model_fit$draws( )[, 1057:1072] ), ncol =  4, byrow =  TRUE)

mS2

fit$model_fit$draws( )[, c( 'S[1,2]','S[1,3]','S[1,4]','S[2,3]','S[2,4]','S[3,4]')]


qgraph::qgraph(mS2, labels = varnames)

## needs to be list of full cormats
qgraph::qgraph.animate( list())
