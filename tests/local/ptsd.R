## Data from: https://osf.io/emwr9/files/
## Howe, Reeves, & Fisher
options(width = 200, pillar.subtle = FALSE)

ptsd <- read.csv( file = "./2022_03_18_MCFA_data.csv" )
head( ptsd )

ptsd_comp <- ptsd#[complete.cases( ptsd ), ]
head( ptsd_comp )

## Select some variables
dat1 <- ptsd_comp[, c('id', 'enthus', 'fear', 'angry', 'happy', 'positive', 'horror')]
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
              v.names = c('enthus', 'fear', 'angry', 'happy', 'positive', 'horror'))

X0
## Replace missings with 50
#X0[is.na(X0 )] <- 50



X0$id <- seq_len(nrow(X0))

XL <- reshape(X0, direction = "long",
              idvar = "id", sep = ".",
              varying = names(X0)[2:length(names(X0))] )

XL
## sort by id, then time
X2 <- XL[order( XL$id,  XL$time), ]
head(X2)

##
X2 <- X2[X2$time <= 70,]
sum(is.na(X2))

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
tsdat <- lapply( seq_len(N), function(x) X2[X2$id == x, 3:(2+nts)])
str(tsdat)

tsdat

devtools::load_all( )

getwd( )
setwd( "../")
tsdat

fit <- dcnet( data = tsdat, J =  N, group =  groupvec, standardize_data = TRUE,
             sampling_algorithm = 'HMC')

summary(fit)


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


fit_const <- dcnet( data = tsdat, J =  N, group =  groupvec, parameterization = 'CCC',
                   standardize_data = TRUE,
                   sampling_algorithm = 'variational')


log_lik_loc <- grep( "log_lik", colnames(fit_const$model_fit$draws( )) )
log_lik_const <- fit_const$model_fit$draws( )[, log_lik_loc]
waic( log_lik = fit_const$model_fit$draws( )[, log_lik_loc] )

loo::waic( log_lik )
m1 <- loo::loo( log_lik )
m2 <- loo::loo( x = log_lik_const)

loo::loo_compare(m1,  m2 )

rts_out_loc <- grep( "rts_out", colnames(fit$model_fit$draws( )) )

rts_out <- fit$model_fit$draws( )[, rts_out_loc]
head(rts_out)
## Format of variable names: rts_out[j, T, nt ]
rts_out[, 18:25]

#saveRDS( fit,  file = './local/ptsd.Rds')
load(file = './local/ptsd.Rds')

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
#'enthus', 'fear', 'angry', 'happy'
c12 <- ggplot(df,  aes(x = time,  y = cor12 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +  scale_y_continuous("PCor(enthus, fear)")
c13 <- ggplot(df,  aes(x = time,  y = cor13 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +  scale_y_continuous("PCor(enthus, angry)")
c14 <- ggplot(df,  aes(x = time,  y = cor14 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) )+  scale_y_continuous("PCor(enthus, happy)")


plout2 <- sapply(3:nts,  function(p) {
  sapply(1:N,  function(f) {
    sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('R[',f,',',x,',2,',p,']')]))
  })
})

df2 <- data.frame( plout2, time = rep(seq(1:tl), N ), id = as.factor(rep(1:N,  each = tl)))
names(df2)[1:2] <- c('cor23',  'cor24' )
head(df2)

c23 <- ggplot(df2,  aes(x = time,  y = cor23 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +  scale_y_continuous("PCor(fear, angry)")
c24 <- ggplot(df2,  aes(x = time,  y = cor24 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +  scale_y_continuous("PCor(fear, happy)")


plout3 <- sapply(4,  function(p) {
  sapply(1:N,  function(f) {
    sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('R[',f,',',x,',3,',p,']')]))
  })
})

df3 <- data.frame( plout3, time = rep(seq(1:tl), N ), id = as.factor(rep(1:N,  each = tl)))
names(df3)[1] <- c('cor34')
head(df3)

c34 <- ggplot(df3,  aes(x = time,  y = cor34 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) )+  scale_y_continuous("PCor(angry, happy)")

nn <- ggplot( ) + theme_void()

ggplot(df3, aes(x = time,  y = cor34 , color = id)) + geom_line( )

library(patchwork )

(c12 | c13 | c14 ) /
( nn | c23 | c24 ) /
  ( nn | nn  | c34 )

ggsave(filename = "/home/philippe/UZH/Kongresse/Kongresse2022/Talk/Figures/pcor.pdf", width = 8, height = 8)



X <- as.matrix( scale( tsdat[[1]], scale = F) )
X

V <- (t(X)%*%X / 69)

D <- matrix(0, ncol = 4, nrow = 4)
diag(D) <- sqrt(diag(V))
Ds <- solve(D)

Ds%*%V%*%Ds

