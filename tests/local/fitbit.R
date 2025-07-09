#####################
## Siwei's fitbit  ##
####################

options(width = 190)
library(data.table)

#fitbit <- data.table(read.csv(file = "~/Dropbox/Projekte/M_GARCH/WORK/-DATASETS/SOURCE/fit_bit_dat.csv")) 
fitbit <- data.table(read.csv(file = "~/SiweisFitbit/fit_bit_dat.csv"))

dim(fitbit)

fitbit_summary <- summary(fitbit)
print(fitbit_summary)

head(fitbit)
str(fitbit)


fitbit$sedentaryHours <- fitbit$sedentaryMinutes/60


## Select variables of interest
variables <- c("active", "totalDistance", "interested", "excited", "attentive", "inspired", "proud", "enthusiastic", "strong", "sedentaryMinutes")
variables <- c("active", "attentive", "sedentaryHours", "totalDistance")


##  Overwrite disinterested with reversing "inerested"

variables <- c("active", "excited", "totalDistance", "interested")
fitcomp <- fitbit[complete.cases(fitbit[, ..variables])]

fitcomp

## Select only individuals with 80+ entries
tmp <- fitcomp[, .N, by = record_id]
tmp
sel <- tmp[N >= 98]$record_id
length(sel)

## Subset selected
dt <- fitcomp[record_id %in% sel]
dt

## Transform the four vaiables to an unbounded space
## dt[, (variables) := lapply(.SD, scale), .SDcols = variables]

logit <- function(x) {
  epsilon <- 1e-1  # Small constant to avoid boundary issues
  yprime <- (x / 100) + epsilon
  yprime <- pmin(pmax(yprime, epsilon), 1 - epsilon)  # Keep yprime within (0, 1)
  logit_y <- log(yprime / (1 - yprime))
  return(logit_y)
}

#dt[, (c("active","excited","disinterested")) := lapply(.SD, logit), .SDcols = c("active","excited","disinterested")]

logit2 <- function(x) {
  epsilon <- 1e-1  # Small constant to avoid boundary issues
  yprime <- (x / 25) + epsilon
  yprime <- pmin(pmax(yprime, epsilon), 1 - epsilon)  # Keep yprime within (0, 1)
  logit_y <- log(yprime / (1 - yprime))
  return(logit_y)
}

#dt[, sedentaryHours := logit2(sedentaryHours)]

#dt[, totalDistance := log(sedentaryMinutes/60 + 1e-1)]

range(dt[, totalDistance])

dt[, ..variables]

## remove data to fre up memory
rm(list = c("fitbit", "fitcomp"))

## Obtain minimal length of time-series (ts)
tsl <- min(dt[, .N, by = record_id]$N)

## Truncate the data table to exactly tsl observations per person
## The .SD (SubsetData) refers to the current group of rows being processed.
dt <- dt[, head(.SD, tsl), by = record_id]

## Re-Create observation index
dt[, time := 1:.N, by = record_id]

# Create a new id variable that is a sequence
dt[, id := .GRP, by = record_id]

## drop 4 individuals with basically no variation
## 23, 28, 31, 34: Didn't seem to wear fitbit all the time
## 10., 17, 20: zero pattern
#dt <- dt[!id %in% c(10, 17, 20, 23, 28, 31, 34)]
#dt <- dt[!id %in% c(11, 40, 46, 48, 50, 51, 52, 54)]

dt <- dt[!id %in% c(11, 12, 33, 36, 40, 45, 46, 48, 50, 51, 52, 54, 55,  9, 24, 32, 36, 53 )]

dt

## Re-Create observation index
dt[, id := .GRP, by = id]
unique(dt$id)

library(ggplot2 )
#ggplot(dt,  aes(x = time,  y = excited, group = factor(id), color = factor(id)) ) +
#  geom_smooth( show.legend = FALSE, se = FALSE) + geom_point(alpha = .2)
head(dt )
unique(dt$id )

plt <- ggplot(dt[id == 2],  aes(x = time,  y = active, group = factor(id), color = factor(id)) ) +
  geom_smooth( show.legend = FALSE, se = FALSE) + geom_point(alpha = .2)
#plt

N <- length(unique(dt$id))
N

variables
plt <- ggplot(dt[id <= N],  aes(x = time,  y = active, group = factor(id), color = factor(id)) ) +
  geom_smooth( show.legend = FALSE, se = FALSE) + geom_point(alpha = .4)
#plt + facet_wrap(~ id)


#ggsave(filename = "~/Nextcloud/Documents/fitbit_raw.pdf", width = 15, height = 15 )
dev.off()


ids <- unique(dt$id)
tl <- tsl
groupvec <-
  rep(seq_len(N),  each = tl)

variables

tsdat <- lapply(ids, function(x) dt[dt$id == x, .SD, .SDcols = variables])
tsdat

## use tsdat to create a similar dataset for use in the Readme and as example data in the package:
## We're going to add jitter to create a smilar copy of the original data

jtr <- function(x) {
  apply(x,  2,  function(slice) {
    jitter( x = unlist(slice), factor = 10)  
  })
}

ema_fitbit <- lapply(seq_len(N),  function(i) {
  jtr(tsdat[[i]])
})

#save(ema_fitbit,  file = "data/ema-fitbit.rda")



getwd()
#setwd( "../")

devtools::load_all(path = "./")

## Model with two S matrices pre/post intervention
fit_init <- dcnet(data = tsdat, J = N, group = groupvec, S_pred = NULL,
                  standardize_data = TRUE,
                  parameterization = "DCCrs",
                  iterations = 30000,
                  sampling_algorithm = "variational",
                  meanstructure = "VAR",
                  chains = 4,
                  init = 0.1)

summary(fit_init)


fit <- dcnet(data = tsdat, J = N, group = groupvec, S_pred = NULL,
             parameterization = "DCCrs",
             iterations = 30000,
             standardize_data = TRUE,
             sampling_algorithm = "variational",
             algorithm = "fullrank",
             #grad_samples = 20,
             #elbo_samples = 200,
             #adapt_iter = 200,
             #eta = 0.005,
             chains = 5,
             init = fit_init$model_fit)

summary(fit)

grep("mu", colnames(fit$model_fit$draws()))
colMeans(fit$model_fit$draws()[, 139346:139350])

grep("phi0_fixed", colnames(fit$model_fit$draws()))
colMeans(fit$model_fit$draws()[, 3:6])

grep("phi0_fixed", colnames(fit$model_fit$draws()))
int <- matrix(colMeans(fit$model_fit$draws()[, 3:6]), ncol = 1)
int

grep('vec_phi_fixed', colnames(fit$model_fit$draws( )) )
phi <- matrix(colMeans(fit$model_fit$draws( )[, 183:198]), ncol = 4)
phi

int + phi %*% t(tsdat[[1]][1,])



#fit <- read_cmdstan_csv( "/tmp/Rtmpv70rbe/DCCMGARCHrandS-202306221623-1-675102.csv", variables = "H")
#fit$draws
                                        #posterior::summarise_draws(fit$draws )

fitC <- dcnet( data = tsdat, J = N, group = groupvec, S_pred = NULL, parameterization = "CCC",
             standardize_data = FALSE, sampling_algorithm = 'variational', threads = 1, init = 0.5, tol_rel_obj = 0.0025)



#summary(fitC)
grep('rtscheck', colnames(fitC$model_fit$draws( )) )
fitC$model_fit$draws( )[, 35372:35375]

grep('phi0_fixed', colnames(fitC$model_fit$draws( )) )
int <- matrix(colMeans(fitC$model_fit$draws( )[, 3:6]), ncol = 1)

grep('vec_phi_fixed', colnames(fitC$model_fit$draws( )) )
phi <- matrix(colMeans(fitC$model_fit$draws( )[, 183:198]), ncol = 4)

int + phi %*% t(tsdat[[1]][1,])


grep('phi0', colnames(fitC$model_fit$draws( )) )
colMeans(fitC$model_fit$draws( )[, grep('phi0', colnames(fitC$model_fit$draws( )) )])



log_lik_r <- grep( "log_lik", colnames(fit$model_fit$draws( )) )
log_lik_r <- fit$model_fit$draws( )[, log_lik_r]
dim(log_lik_r )

r_eff_r <- loo::relative_eff(exp(log_lik_r ),  chain_id = rep(1,  1000 ) )

log_lik_c <- grep( "log_lik", colnames(fitC$model_fit$draws( )) )
log_lik_c <- fitC$model_fit$draws( )[, log_lik_c]
dim(log_lik_c )
r_eff_c <- loo::relative_eff(exp(log_lik_c ),  chain_id = rep(1,  1000 ) )


fr <- loo::loo(log_lik_r, r_eff = r_eff_r)
fc <- loo::loo(log_lik_c, r_eff = r_eff_c )
fr
fc
loo::loo_compare(fr,  fc )
## random efx model is preferred

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
abbreviate(variables)
varnames <- abbreviate(variables, minlength = 2)
varnames
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

ggsave(filename = "~/Dropbox/Public/singlecorr.pdf",  width = 6,  height = 3,  plot = c34 )

nn <- ggplot( ) + theme_void()

library(patchwork )

patched <- (c12 | c13 | c14 ) /
( nn | c23 | c24 ) /
( nn | nn  | c34 )

patched <- (c12 | c13 | c14 ) /
( c34 | c23 | c24 )
patched

ggsave(filename = "~/Nextcloud/Documents/pcor.pdf", width = 9, height = 4.5, plot = patched)


### rtsout holds yrep
##Yrep => rts_out[N,ts_length,timesries]


fit$model_fit$draws( )[,'rts_out[34,68,4]'][1:5,]

## Obtian random rows with 99 observations across participants (row contains the values for the time series 99; there are a total of 99 x 1000 samples. Column contains subject)

variable <- 4
variables[variable]

yrep1 <- sapply(seq_len(N), function(person) {
  sapply(seq_len(tl), function(x) {
    fit$model_fit$draws( )[, paste0('rts_out[',person,',',x,',',variable,']')]
  })
})

head(yrep1)

tsl
yrep1ext <- cbind(yrep1, 1:tsl, rep(1:1000,  each = tsl ))

## randomly select 100 samples from yrep1ext
sel <- sample(1:1000,  100 )
## then paste these 100 samples into one row, so that we have a matirx with 100 rows representing the samples and the columns concatenate all 34 indivduals time their 99 time points

ysamp <- list()
for(i in 1:100) {
  smp <- sel[i]
  ysamp[[i]] <- c(yrep1ext[yrep1ext[,N+2] == smp, c(1:N)])
}

## Convert the list into a matrix with 100 rows
ysamp_matrix <- do.call(rbind, ysamp)

## Check the dimensions of the matrix
dim(ysamp_matrix)

variables
plots <- list( )
for(i in 1:N ) {
  p <- bayesplot::ppc_dens_overlay(y = dt[dt$id == i, get(variables[variable])],  yrep = ysamp_matrix[,(1:tsl)+(i-1)*tsl], trim =  TRUE )
  plots[[i]] <- p
}

ppc_wrap <- wrap_plots(plots,  ncol = 6 )
ppc_wrap

ggsave(filename =paste0("~/Nextcloud/Documents/ppc_wrap_", variables[variable],".pdf"), width = 15, height = 15, plot = ppc_wrap)

#plots[[2]]
#plots[[13]]
#plots[[20]]
#plots[[33]]

## PPC across N
dim(ysamp_matrix )
length(dt[, get(variables[variable])])
ppc <- bayesplot::ppc_dens_overlay(y = dt[, get(variables[variable])],  yrep = ysamp_matrix , trim = TRUE)

ppc
ppcecd <- bayesplot::ppc_ecdf_overlay( y = dt[, get(variables[variable])],  yrep = ysamp_matrix,
                                      scale_x_continuous(limits = c(0, 100)))
ppcecd

ggsave(filename = paste0("~/Nextcloud/Documents/ppc_", variables[variable], ".pdf"), width = 9, height = 4.5, plot = ppc)
ggsave(filename = paste0("~/Nextcloud/Documents/ppcecd_", variables[variable], ".pdf"), width = 9, height = 4.5, plot = ppcecd)

## Individual predictions
sc <- list()
for(i in 1:N) {
  person <- i
  ## extract lower, median, upper quantile and scal back to original scale
  yrep4 <- sapply(seq_len(tl), function(x) {
    quantile(fit$model_fit$draws( )[,paste0('rts_out[',person,',',x,',',variable,']')], c(.025, .5, .975))
  })
  yrep4 <- data.frame(yrep = t(yrep4))
  yrep4$time <- seq_len(tl )
  yrep4$obs <- unlist(tsdat[[person]][,..variable])
  sc[[i]] <- ggplot(yrep4,  aes(x = time, y = yrep.50.) ) + geom_line( ) + geom_ribbon(aes(ymin =  yrep.2.5.,  ymax = yrep.97.5.),  alpha = .2) + geom_line( aes(y = obs ), color = 'red')+ylab(varnames[variable])
}

ppcind <- wrap_plots(sc,  ncol = 6 )
ppcind

ggsave(filename = paste0("~/Nextcloud/Documents/ppcind_", variables[variable], ".pdf"), width = 15, height = 9, plot = ppcind)


#ggsave(filename = "~/Dropbox/Public/sancheck4.pdf",  plot = sc ,  width = 6, height = 4)


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

sc <- ggplot(yrep4,  aes(x = time, y = yrep.50.) ) + geom_line( ) + geom_ribbon(aes(ymin =  yrep.2.5.,  ymax = yrep.97.5.),  alpha = .2) + geom_line( aes(y = obs ), color = 'red')+ylab(varnames[1])
sc
ggsave(filename = "~/Dropbox/Public/sancheck1.pdf",  plot = sc ,  width = 6, height = 4)


## Sanity check to see if stan_model transform to original data
#yrep4$standat <- fit$RTS_full[[person]][,4] #*fit$grand_sd[1] + fit$grand_mean[1]
#ggplot(yrep4,  aes(x = time, y = standat) ) + geom_point( ) + geom_ribbon(aes(ymin =  yrep.2.5.,  ymax = yrep.97.5.),  alpha = .2) + geom_line( aes(y = obs ), color = 'red')

library(qgraph )

minsf <- min((grep('Sfixed', colnames(fit$model_fit$draws( )) )))
colMeans( fit$model_fit$draws( )[, minsf:(minsf+15)] )

mS <- matrix(colMeans( fit$model_fit$draws( )[, minsf:(minsf+15)] ), ncol =  4, byrow =  TRUE)
fit$model_fit$draws( )[, c( 'Sfixed[1,2]','Sfixed[1,3]','Sfixed[1,4]','Sfixed[2,3]','Sfixed[2,4]','Sfixed[3,4]')]

pdf(file = "~/Nextcloud/Documents/qgraph.pdf")
gR <- qgraph::qgraph(mS, labels = varnames)
dev.off()


## Plot individual Random offects of S
namemat <- matrix(apply(expand.grid( varnames,  varnames ), 1, FUN = function(x) paste0(x[1],"-", x[2])), ncol = 4)
namevech <- namemat[lower.tri(namemat)]
indR <- data.frame(correlation = NA, namevech = NA, subject = NA)
indR

for(i in seq_len(N) ) {
  pos <-  (grep(paste0('S\\[',i,','), colnames(fit$model_fit$draws( )) ))
  #smat <- matrix(colMeans( fit$model_fit$draws( )[, pos ]), ncol = 4, byrow = TRUE)
  smat <- matrix(apply(fit$model_fit$draws( )[, pos ], 2, median), ncol = 4, byrow = TRUE)
  indR <- rbind(indR, data.frame(correlation = smat[lower.tri(smat)], namevech, subject = i))  
}
indR <-  data.frame( indR[-1,] )

mean_data <- data.frame( namevech, mean_value = mS[lower.tri(mS)])
mean_data

rind <- ggplot( indR,  aes(x = correlation, y = subject, color = as.factor(subject)) ) +
  geom_vline(xintercept = 0)+
  geom_vline(data = mean_data, aes(xintercept = mean_value),
             color = "gray")+
  geom_point(show.legend = FALSE ) +
  xlim(-0.1,  0.6)+
  facet_wrap(. ~ namevech, nrow = 2)+
  xlab("Unconditional Correlation" ) +
  ylab("Individual")
rind

ggsave(filename = "~/Nextcloud/Documents/Rind.pdf", width = 8, height = 8, plot = rind)


# Contemporaneous graphs:
## contempraneous model for VAR
minsf <- min(grep('rescov', colnames(fitC$model_fit$draws( )) ))
colMeans( fitC$model_fit$draws( )[, minsf:(minsf+15)] )

mT <- matrix(colMeans( fitC$model_fit$draws( )[, minsf:(minsf+15)] ), ncol =  4, byrow =  TRUE)
mT

## Convert Covariance matrix to partial correlation
theta <- solve(mT )
Delta <- diag( 1 / sqrt( diag(theta) ) )
## Partial Correlation
pc <- - Delta %*% theta %*% Delta + 2*diag(dim(mT)[1])
pc


#pdf(file = "~/Nextcloud/Documents/qgraph_var_only.pdf")
gC <- qgraph::qgraph(round(pc, 5), labels = varnames); mtext("    Panel B", side = 3, adj = 0, line = 1.5, cex = 1.2)
                                        #dev.off()



pdf(file = "~/Nextcloud/Documents/qgraph_RC.pdf", width = 9, height = 5)
layout(matrix(c(1,2), 1,2) )
layout <- qgraph(mS, labels = varnames); mtext("    Panel A", side = 3, adj = 0, line = 1.5, cex = 1.2)
layout <- qgraph(round(pc, 5), labels = varnames); mtext("    Panel B", side = 3, adj = 0, line = 1.5, cex = 1.2)
dev.off()

# Temporal Graphs:
## Temporal graph for VAR-DCC
misfc <- min( (grep('phi_fixed', colnames(fit$model_fit$draws( )) )) )
misfc
colMeans( fit$model_fit$draws( )[, misfc:(misfc)] )

temporal_m <- matrix(colMeans( fit$model_fit$draws( )[, misfc:(misfc+15)] ), ncol =  4, byrow =  FALSE)
temporal_m

pdf(file = "~/Nextcloud/Documents/temporal.pdf")
qgraph::qgraph(temporal_m, labels = varnames)
dev.off()

## Temporal graph for VAR only
(grep('phi_fixed', colnames(fitC$model_fit$draws( )) ))

colMeans( fitC$model_fit$draws( )[, misfc:(misfc+15)] )

temporal_mc <- matrix(colMeans( fitC$model_fit$draws( )[, misfc:(misfc+15)] ), ncol =  4, byrow =  FALSE)
temporal_mc

pdf(file = "~/Nextcloud/Documents/temporal_var.pdf")
qgraph::qgraph(temporal_mc, labels = varnames)
dev.off()



## Precompile models:
library(cmdstanr )
stan_path <- "~/Dropbox/Git/dcnet/inst/stan/"
dcc_file <- file.path(stan_path,  "DCCMGARCHrandQ.stan")
cc_file <- file.path(stan_path, "VAR.stan")
dccr_file <- file.path(stan_path, "DCCMGARCHrandS.stan")

cmdstan_model(dcc_file, dir = stan_path, force_recompile = TRUE)$save_model_file(file = file.path(stan_path, "DCCMGARCHrandQ"))
cmdstan_model(cc_file, dir = stan_path, force_recompile = TRUE)$save_model_file(file = file.path(stan_path, "VAR"))
cmdstan_model(dccr_file, dir = stan_path, force_recompile = TRUE)$save_model_file(file = file.path(stan_path, "DCCMGARCHrandS"))
