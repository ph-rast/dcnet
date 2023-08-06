#####################
## Siwei's fitbit  ##
####################
options(width = 210)
library(data.table )

fitbit <- data.table( read.csv(file = '~/Dropbox/Projekte/M_GARCH/WORK/-DATASETS/SOURCE/fit_bit_dat.csv') )
dim(fitbit )

head(fitbit)
str(fitbit )

## Select variables of interest
variables <- c("excited", "disinterested", "active", "totalDistance")
fitcomp <- fitbit[complete.cases( fitbit[, ..variables] )]

## Select only individuals with 80+ entries
tmp <- fitcomp[, .N, by = record_id]
tmp
sel <- tmp[N >= 95]$record_id
length(sel)

## Subset selected
dt <- fitcomp[record_id %in% sel]

## remove data to fre up memory
rm( list = c("fitbit", "fitcomp"))

## Obtain minimal length of time-series (ts)
tsl <- min(dt[, .N, by = record_id]$N)

## Truncate the data table to exactly tsl observations per person
## The .SD refers to the current group of rows being processed.
dt <- dt[, head(.SD, tsl), by = record_id]

## Re-Create observation index
dt[, time := 1:.N, by = record_id]

# Create a new id variable that is a sequence
dt[, id := .GRP, by = record_id]


library(ggplot2 )
ggplot(dt,  aes(x = time,  y = excited, group = factor(id), color = factor(id)) ) +
  geom_smooth( show.legend = FALSE, se = FALSE) + geom_point(alpha = .2)

dev.off( )

N <- length(unique(dt$id ) )
N

ids <- unique(dt$id)
tl <- tsl
groupvec <- rep(seq_len(N),  each = tl)

variables

tsdat <- lapply(ids, function(x) dt[dt$id == x, .SD, .SDcols = variables] )
tsdat


getwd( )
setwd( "./tests/")

devtools::load_all( )

## Model with two S matrices pre/post intervention
fit <- dcnet( data = tsdat, J = N, group = groupvec, S_pred = NULL, parameterization = "DCCr",
             standardize_data = FALSE, sampling_algorithm = 'variational', threads = 1, init = 1)


#fit <- read_cmdstan_csv( "/tmp/Rtmpv70rbe/DCCMGARCHrandS-202306221623-1-675102.csv", variables = "H")
#fit$draws
#posterior::summarise_draws(fit$draws )



fit$model_fit$draws(variables = "pcor[1,1,1,1]" )

## Plots
ind <- seq_len(N)
nts <- 4

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
varnames <- variables
#'enthus', 'fear', 'angry', 'happy'
c12 <-
  ggplot(df,  aes(x = time,  y = cor12 , color = id)) + geom_line(show.legend = FALSE )+ coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous( paste0("PCor(", varnames[1], ", ", varnames[2], ")") ) #+ geom_vline( xintercept = 36)
c13 <-
  ggplot(df,  aes(x = time,  y = cor13 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous( paste0("PCor(", varnames[1], ", ", varnames[3], ")") )#+ geom_vline( xintercept = 36)
c14 <- ggplot(df,  aes(x = time,  y = cor14 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous( paste0("PCor(", varnames[1], ", ", varnames[4], ")") )#+ geom_vline( xintercept = 36)



df2 <- data.frame( plout2, time = rep(seq(1:tl), N ), id = as.factor(rep(ind,  each = tl)))
names(df2)[1:2] <- c('cor23',  'cor24' )
head(df2)

c23 <- ggplot(df2,  aes(x = time,  y = cor23 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous( paste0("PCor(", varnames[2], ", ", varnames[3], ")") )#+ geom_vline( xintercept = 36)
c24 <- ggplot(df2,  aes(x = time,  y = cor24 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous( paste0("PCor(", varnames[2], ", ", varnames[4], ")") )#+ geom_vline( xintercept = 36)



df3 <- data.frame( plout3, time = rep(seq(1:tl), N ), id = as.factor(rep(ind,  each = tl)))
names(df3)[1] <- c('cor34')
head(df3)

c34 <- ggplot(df3,  aes(x = time,  y = cor34 , color = id)) + geom_line(show.legend = FALSE ) + coord_cartesian(ylim = c(-1, 1 ) ) +
  scale_y_continuous(paste0("PCor(", varnames[3], ", ", varnames[4], ")"))#+ geom_vline( xintercept = 36)

nn <- ggplot( ) + theme_void()

library(patchwork )

(c12 | c13 | c14 ) /
( nn | c23 | c24 ) /
( nn | nn  | c34 )


### Yrep => rts_out[N,ts_length,timesries]
tsl
N
fit$model_fit$draws( )[,'rts_out[46,68,4]']

person <- 25
## extract lower, median, upper quantile and scal back to original scale
yrep4 <- sapply(seq_len(tl), function(x) {
  quantile(fit$model_fit$draws( )[, paste0('rts_out[',person,',',x,',1]')], c(.025, .5, .975))
  })

yrep4

#yrep4 <- yrep4*fit$grand_sd[1] + fit$grand_mean[1]

yrep4 <- data.frame(yrep = t(yrep4))


yrep4$time <- seq_len(tl )
names(yrep4)
yrep4$obs <- unlist( tsdat[[person]][,1] )

tsdat[[1]][,1]
head(yrep4 )

ggplot(yrep4,  aes(x = time, y = yrep.50.) ) + geom_line( ) + geom_ribbon(aes(ymin =  yrep.2.5.,  ymax = yrep.97.5.),  alpha = .2) + geom_line( aes(y = obs ), color = 'red')


## Sanity check to see if stan_model transform to original data
yrep4$standat <- fit$RTS_full[[person]][,1] #*fit$grand_sd[1] + fit$grand_mean[1]
ggplot(yrep4,  aes(x = time, y = standat) ) + geom_point( ) + geom_ribbon(aes(ymin =  yrep.2.5.,  ymax = yrep.97.5.),  alpha = .2) + geom_line( aes(y = obs ), color = 'red')

