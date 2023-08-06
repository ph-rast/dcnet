#####################
## Siwei's fitbit  ##
####################

library(data.table )

fitbit <- data.table( read.csv(file = '~/Dropbox/Projekte/M_GARCH/WORK/-DATASETS/SOURCE/fit_bit_dat.csv') )
dim(fitbit )

head(fitbit)

## Select variables of interest
variables <- c("excited", "upset", "enthusiastic", "stressed")
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
setwd( "../")

devtools::load_all( )

## Model with two S matrices pre/post intervention
fit <- dcnet( data = tsdat, J = N, group = groupvec, S_pred = NULL, parameterization = "DCCr",
             standardize_data = TRUE, sampling_algorithm = 'variational', threads = 1, init = 1)


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
