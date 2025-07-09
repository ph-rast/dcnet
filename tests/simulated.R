##library(dcnet )
# we recommend running this is a fresh R session or restarting your current session
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#remotes::install_github("stan-dev/cmdstanr")
#cmdstanr::install_cmdstan( )

set.seed(234)

ivv <- function(V, nt) {
    z <- diag(0, nt)
    L <- diag(0, nt)
    index <- 0
    for (j in 1:nt) {
        for (i in 1:nt) {
            if (i > j) {
                index <- index + 1
                z[i, j] <- V[index]
            }
            if (i < j) {
                z[i, j] <- 0
            }
        }
    }
    for (i in 1:nt) {
        for (j in 1:nt) {
            if (i < j) {
                L[i, j] <- 0
            }
            if (i == j) {
                if (i == 1) {
                    L[i, j] <- 1
                }
                if (i > 1) {
                    L[i, j] <- sqrt(1 - sum(L[i, 1:j]^2))
                }
            }
            if (i > j) {
                L[i, j] <- z[i, j] * sqrt(1 - sum(L[i, 1:j]^2))
            }
        }
    }
    L
}


L <- ivv(tanh(c(0.3, -0.5, 1.1, 0, 0.3, -.2)), 4)
L

L %*% t(L)

devtools::load_all()
options(width = 200, crayon.enabled = FALSE)

S <- list()
for(j in 1:3) {
  S[[j]] <- dcnet:::.ranS(4)
}
S[[1]]
array(S, c(4,4,3) )


## Create data:
N <- 40
tl <- 50
nts <- 4
simdat <- .simuDCC(tslength = tl,  N = N,  n_ts = nts,
                   phi0_sd = 0.5,  ## random effects
                   log_c_fixed = rep(1, nts),
                   log_c_r_sd = 0.5,
                   a_h_fixed = rep(-1.5, nts),
                   a_h_r_sd = 0.5,
                   b_h_fixed = rep(-0.5, nts),  ## On logit scale
                   b_h_r_sd = 0.5,
                   l_a_q_fixed = -1.5,  ## on logit scale
                   l_b_q_fixed = -0.5,   ## on logit scale
                   l_a_q_r_sd = 0.5,
                   l_b_q_r_sd = 0.5,
                   phi0_fixed =  rep(0, nts),
                   ranS_sd = 0.1, ## random effects on atanh scale
                   phi_mu = 0, ## populate phi
                   phi_sd = 0.1, ## create variation in the intercepts of the time series
                   phi_ranef_sd = 0.1, ## ranef of the parameter matrix
                   stationarity_phi = TRUE)


rtsgen <- lapply(seq(dim(simdat[[1]])[3]), function(x) t( simdat[[1]][,,x] ))
groupvec <- rep(c(1:N),  each = tl)

atanh(simdat$Fixed_S[lower.tri(simdat$Fixed_S)])
simdat$fixed_S_atanh

rtsgen[[1]][,1]

## Fixed Corr
csimdat$S
simdat$Fixed_S
simdat$phi0_fixed
simdat$fixed_phi


person <- 2
plot(rtsgen[[person]][, 1], type = "l")#, ylim = c(-6, 6))
lines(rtsgen[[person]][, 2], lty = 2)
lines(rtsgen[[person]][, 3], lty = 3, col = "red")
#lines(rtsgen[[person]][, 4], lty = 4, col = "blue")

dev.off()
simdat$fixed_phi

init_fun <- function() {
  list(
    c_h        = rep(0.1, nts),          # vector[nt]
    a_h_raw    = rep(0.05, nts),         # if you use raw scale
    b_h_raw    = rep(0.90, nts),
    phi0_fixed = rep(0, nts),
    vec_phi_fixed = rep(0, nts * nts),
    # add any other required parameters here...
    phi_stdnorm = matrix(0, N, nts * nts)    # zeros for subject draws
  )
}

devtools::load_all()

system.time({

    fit0 <- dcnet(
        data = rtsgen, parameterization = "CCC", J =  N,
        group = groupvec, standardize_data = FALSE,
        init = 0,
        meanstructure = "VAR",
        iterations = 30000,
        chains = 4,
        threads = 5, tol_rel_obj =  0.005, ## 8 threads: 188 mins /
        sampling_algorithm = "pathfinder", # "laplace",
        grainsize = 10)

    
    fit0 <- dcnet(
        data = rtsgen, parameterization = "DCCrs", J =  N,
        group = groupvec, standardize_data = FALSE,
        init = 0,
        meanstructure = "VAR",
        iterations = 30000,
        chains = 4,
        threads = 4, tol_rel_obj =  0.005, ## 8 threads: 188 mins /
        sampling_algorithm = "variational",#"pathfinder", # "laplace",
        grainsize = 3)

    fit <- dcnet(
        data = rtsgen, parameterization = "DCCrs", J = N,
        group = groupvec, standardize_data = FALSE,
        init = fit0$model_fit,
        meanstructure = "VAR",
        iterations = 1000,
        sampling_algorithm = "hmc",
#        algorithm = "fullrank", ## fullrank should be less biased
        #grad_samples = 1,
        #elbo_samples = 150,
        #eta = 0.25,
        #adapt_iter = 500,
        #grainsize = 3,
        #threads_per_chain = 4
        chains = 4
    ) ## grain is subject. Make it max 3 per chunk, then grainsize * threads =(approx) subjects
    ## 46 mins to beat (gs=30) / gs=1, 60 mins / gs = 2, 118 / gs=3, 93 mins / gs=4,  / gs=6,   
    ##num_threads =  8)

})



fit$model_fit$output()


fit

######################################
## Simulate 
######################################

## list containing all replicaion datasets:
replication_data <- list()

## Create data:
simulate_data <- function(N = 15, tl = 50, nts = 3) {
    simdat <- .simuDCC(
        tslength = tl, N = N, n_ts = nts,
        phi0_sd = 0.5, ## random effects
        log_c_fixed = rep(1, nts),
        log_c_r_sd = 0.25,
        a_h_fixed = rep(-1.5, nts),
        a_h_r_sd = 0.5,
        b_h_fixed = rep(-0.5, nts), ## On logit scale
        b_h_r_sd = 0.5,
        l_a_q_fixed = -1.5, ## on logit scale
        l_b_q_fixed = -0.5, ## on logit scale
        l_a_q_r_sd = 0.5,
        l_b_q_r_sd = 0.5,
        phi0_fixed = rep(0, nts),
        ranS_sd = 0.25, ## random effects on atanh scale
        phi_mu = 0, ## populate phi
        phi_sd = 0, ## create variation in the intercepts of the time series
        phi_ranef_sd = 0.1, ## ranef of the parameter matrix
        stationarity_phi = FALSE
    )
    rtsgen <- lapply(seq(dim(simdat[[1]])[3]), function(x) t(simdat[[1]][, , x]))
    groupvec <- rep(c(1:N), each = tl)
    return(list(rtsgen, groupvec, simdat, N = N))
}


replication_data <- simulate_data(N = 10)
replication_data[[1]]
replication_data[[2]]
replication_data[[3]]
replication_data$N

## Function to prevent loop to stop when stan model encounters errors
## Try fitting model 3 times before moving on

## V2:
safe_sample <- function(s, max_retries = 3, replication_data) {
    replication_data <- replication_data

    ## First run meanfield algo to obtain init values:
    fit_init <- dcnet(
        data = replication_data[[1]], parameterization = "DCCrs", J = replication_data$N,
        group = replication_data[[2]], standardize_data = FALSE,
        init = 0,
        meanstructure = "VAR",
        iterations = 50000,
        eta = 0.05,
        tol_rel_obj =  0.005,
        sampling_algorithm = "variational")
    
    ## Helper function to check if model_fit is broken
    is_model_fit_broken <- function(fit) {
        tryCatch({
            if (!inherits(fit$model_fit, "CmdStanFit")) return(TRUE)
            fit$model_fit$metadata()  # harmless method, fails if object is invalid
            return(FALSE)
        }, error = function(e) {
            return(TRUE)
        })
    }

    ## Fit model with fullrank
    for (attempt in seq_len(max_retries)) {
        message(sprintf("Fitting attempt %d for index %d", attempt, s))

        fit <- tryCatch({
            fit_init
            ## dcnet(
            ##    data = replication_data[[1]],
            ##    parameterization = 'DCCrs',
            ##    J = replication_data$N,
            ##    meanstructure = "VAR",
            ##    group = replication_data[[2]], ## groupvec
            ##    standardize_data = FALSE,
            ##    init = fit_init$model_fit,
            ##    threads = 4,
            ##    iterations = 50000,
            ##    #calgorithm = "fullrank",
            ##    sampling_algorithm = 'variational',
            ##    tol_rel_obj =  0.005,
            ##    eta = 0.05,
            ##    adapt_iter = 200,
            ##    chains = 4,
            ##    grainsize = 3
            ## )
        }, error = function(e) {
            message(sprintf("Hard failure on attempt %d: %s", attempt, e$message))
            return(NULL)
        })

        if (is.null(fit)) next

        ## Check if model_fit is broken
        if (is_model_fit_broken(fit)) {
            message("model_fit is invalid or unusable — retrying...")
            next
        }

        ## Try to extract model output
        output_lines <- tryCatch({
            out <- fit$model_fit$output()
            if (is.character(out)) out else as.character(unlist(out))
        }, error = function(e) {
            message("Could not retrieve model output — assuming failure.")
            return("ERROR_FETCHING_OUTPUT")
        })

        ## Patterns that indicate failure or non-convergence
        failure_patterns <- c(
            "All proposed step-sizes failed",
            "algorithm may not have converged",
            "Exception:.*not positive definite",
            "Exception:",
            "Fitting failed. Unable to print"
        )

        failed <- any(sapply(failure_patterns, function(pat) {
            any(grepl(pat, output_lines, ignore.case = TRUE, perl = TRUE))
        }))

        if (failed) {
            message("Detected Stan warning or internal failure — retrying...")
            next
        }

        ## Everything looks okay
        return(fit)
    }
    message(sprintf("All %d attempts failed for index %d.", max_retries, s))
    return(NULL)
}

## END V2

## Init all objects:
phi0_fixed_cov <- phi0_sd_cov <-phi_fixed_covr <-phi_ranef_covr <-log_c_fixed_covr <-log_c_rand_covr <-log_a_fixed_covr <-a_h_r_covr <-log_b_fixed_covr <-b_h_r_covr <-l_a_q_fixed_covr <-l_a_q_r_covr <-l_b_q_fixed_covr <-l_b_q_r_covr <-fix_S <- ran_S <- list()



## Function to compute coverage probability ranging from L to U in the original population distribution
overlap <- function(population, L, U) {
    coverage <- as.numeric(population > L & population < U)
    return(coverage)
}

rmse <- function(model, population) {
    rmse <- sqrt(mean((model - population)^2))
    return(rmse)
}

bias <- function(model, population) {
    bias <- mean(model - population)
    return(bias)
}

#variables

sbc <- function(model, population, column) {
  ## Algorithm 1 from: Validating Bayesian Inference Algorithms with Simulation-Based Calibration (2020)
  bin <- list()
  binl <- 0
  ## Assuming 1000 draws in the fit objects (standard)
  ## This creates 20 bins 
  for(start in seq(1, 951, by=50)) {
    binl <- binl+1
    end <- start + 49
    current_sequence <- start:end
    bin[binl] <- sum(model[start:end,column] < population[start:end,column])
  }
  return(bin)
}


#sbc(fit_r$model_fit$draws(variables = 'phi0_fixed'), fit$model_fit$draws(variables = 'phi0_fixed'),  column = 1 )


variables_m <- c('phi0_fixed', 'phi0_tau', 'vec_phi_fixed', 'phi_tau', 'c_h_fixed', 
               'c_h_tau', 'a_h_fixed', 'a_h_tau', 'b_h_fixed', 'b_h_tau', 
               'l_a_q', 'l_a_q_sigma', 'l_b_q', 'l_b_q_sigma', 'S_vec_fixed', 'S_vec_tau')

var <- c(
    "phi0_fixed", "phi0_sd", "fixed_phi", "phi_ranef_sd", "log_c_fixed",
    "log_c_r_sd", "a_h_fixed", "a_h_r_sd", "b_h_fixed", "b_h_r_sd",
    "l_a_q_fixed", "l_a_q_r_sd", "l_b_q_fixed", "l_b_q_r_sd", "fixed_S_atanh", "ranS_sd"
)

# Initialize empty lists to store the results for each variable
cov_list <- list()
rmse_list <- list()
bias_list <- list()
bins_list <- list()


for (s in 1:100) {

    replication_data <- simulate_data(N = 30)
    fit_r <- safe_sample(s, replication_data = replication_data)

    if (is.null(fit_r)) {
        next
    }

    for (p in 1:length(variables)) {
        SF <- fit_r$model_fit$summary(
            variables = variables_m[p], "mean",
            extra_quantiles = ~ posterior::quantile2(., probs = c(0.025, 0.975))
        )
        cov_list[[variables_m[p]]][[s]] <-
            overlap(
                unlist(replication_data[[3]][var[p]]),
                L = SF$q2.5, U = SF$q97.5
            )

        rmse_list[[variables_m[p]]][[s]] <-
            sapply(seq_len(nrow(SF)), function(i) {
                rmse(
                    fit_r$model_fit$draws(variables = variables_m[p])[, i],
                    unlist(replication_data[[3]][var[p]])
                )
            })
        bias_list[[variables_m[p]]][[s]] <-
            sapply(seq_len(nrow(SF)), function(i) {
                bias(
                    fit_r$model_fit$draws(variables = variables_m[p])[, i],
                    unlist(replication_data[[3]][var[p]])
                )
            })
    }
    gc()
    s
}

## Try reducing the convergence criterion to 0.005
cov_list
rmse_list
bias_list
bins_list

cov_list


var_averages <- list()
for (i in 1:length(cov_list)) {
    nested_average <- c()                                             
    for (j in 1:length(cov_list[[i]])) {
        if (is.null(cov_list[[i]][[j]])) cov_list[[i]][[j]] <- NA
        nested_average <- c(nested_average, mean(cov_list[[i]][[j]], na.rm = TRUE))
        var_averages[[i]] <- nested_average
    }
}
names(var_averages ) <- names(cov_list )
var_averages


var_av <- sapply(var_averages, function(x ) mean(x, na.rm = TRUE))
var_av


n30tl50 <- round(var_av, 2)
n30tl50
##n50tl50 <- round(var_av, 2)
n50tl50
## n50tl75 <- round(var_av, 2)
n50tl75
## n75tl75 <- round(var_av, 2)
n75tl75
## n75tl100 <- round(var_av, 2)
n75tl100
## n100tl100 <- round(var_av, 2)
n100tl100



rbind(n30tl50, n50tl50, n50tl75, n75tl75, n75tl100)

write.csv( rbind(n30tl50, n50tl50, n50tl75, n75tl75, n75tl100), file = "VB_SimulationResults.csv")

bin <- list()
binl <- 0
## Assuming 1000 draws in the fit objects (standard) 
for(start in seq(1, 951, by=50)) {
  binl <- binl+1
  end <- start + 49
  current_sequence <- start:end
  bin[binl] <- sum(fit_r$model_fit$draws( variables = variables[1] )[start:end,col] < fit$model_fit$draws(variables = variables[1] )[start:end,col])
}
bin


unlist(bin)
chisq.test(unlist(bin))
hist(unlist(bin), breaks = length(unlist(bin) ), xlim = c(0, 50 ))

dev.off( )



######################################################################################
######################################################################################
## Compare values to R simu
###





svf <- posterior::as_draws_matrix( fit$model_fit$draws(variables = 'S_vec_fixed' ) )
svf
dt <- data.table(fit$model_fit$summary(variables = c('phi0_fixed', 'S_vec_fixed', 'vec_phi_fixed')))[,c(1,2,6,7)]
dt

## Coverage: phi0_fixed
SF <- fit$model_fit$summary(variables = 'phi0_fixed', "mean",
                      extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
SF

phi0_fixed <- simdat$phi0_fixed
phi0_fixed
phi0_fixed_cov <- phi0_fixed >= SF[,3] & phi0_fixed <= SF[,4]
print(phi0_fixed_cov)

## Coverage: phi0_sd
SF <- fit$model_fit$summary(variables = 'phi0_tau',
                      "mean",
                      extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
SF
phi0_sd <- simdat$phi0_sd
phi0_sd
phi0_sd_cov <- phi0_sd >= SF[,3] & phi0_sd <= SF[,4]
phi0_sd_cov

## Phi is transformed in simulation script to maintain stationarity: This might not
## capture the actual coverage -- needs rethinking
## Coverage: phi_fixed
SF <- fit$model_fit$summary(variables = 'vec_phi_fixed',
                      "mean",
                      extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
SF
simdat$fixed_phi
phi_fixed_covr <- c(simdat$fixed_phi) >= SF[,3] & c(simdat$fixed_phi) <= SF[,4]
print(phi_fixed_covr)
mean(phi_fixed_covr)


## Coverage: phi_ranef_sd
SF <- fit$model_fit$summary(variables = 'phi_tau',
                      "mean",
                      extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
SF
phi_ranef_sd = simdat$phi_ranef_sd
phi_ranef_sd
phi_ranef_covr <-  phi_ranef_sd >=  SF[,3] & phi_ranef_sd <= SF[,4]
print(phi_ranef_covr)

## Covr: Log C fixed: log_c_fixed = rep(0, n_ts), ## on log scale
SF <- fit$model_fit$summary(variables = 'c_h_fixed',
                      "mean",
                      extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
SF
log_c_fixed = simdat$log_c_fixed
log_c_fixed
log_c_fixed_covr <- log_c_fixed >= SF[,3] & log_c_fixed <= SF[,4]
print(log_c_fixed_covr)

## Covr: random effect on C
SF <- fit$model_fit$summary(variables = 'c_h_tau',
                      "mean",
                      extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
SF
log_c_r_sd = simdat$log_c_r_sd
log_c_r_sd
log_c_rand_covr <- log_c_r_sd >= SF[,3] & log_c_r_sd <= SF[,4]
print(log_c_rand_covr)

## Covr: Log A fixed: log_a_fixed = rep(0, n_ts), ## on log scale
SF <- fit$model_fit$summary(variables = 'a_h_fixed',
                      "mean",
                      extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
SF
a_h_fixed = simdat$a_h_fixed
a_h_fixed
log_a_fixed_covr <- a_h_fixed >= SF[,3] & a_h_fixed <= SF[,4]
print(log_a_fixed_covr)

## Covr: random effect on A
SF <- fit$model_fit$summary(variables = 'a_h_tau',
                      "mean",
                      extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
SF
a_h_r_sd = simdat$a_h_r_sd
a_h_r_sd
a_h_r_covr <- a_h_r_sd >= SF[,3] & a_h_r_sd <= SF[,4]
print(a_h_r_covr)

## Covr: Log B fixed: log_b_fixed 
SF <- fit$model_fit$summary(variables = 'b_h_fixed',
                      "mean",
                      extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
SF
b_h_fixed = simdat$b_h_fixed
b_h_fixed
log_b_fixed_covr <- b_h_fixed >= SF[,3] & b_h_fixed <= SF[,4]
print(log_b_fixed_covr)

## Covr: random effect on b
SF <- fit$model_fit$summary(variables = 'b_h_tau',
                      "mean",
                      extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
SF
b_h_r_sd = simdat$b_h_r_sd
b_h_r_sd
b_h_r_covr <- b_h_r_sd >= SF[,3] & b_h_r_sd <= SF[,4]
print(b_h_r_covr)


## Covr: l_a_q 
SF <- fit$model_fit$summary(variables = 'l_a_q',
                      "mean",
                      extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
SF
l_a_q_fixed = simdat$l_a_q_fixed
l_a_q_fixed
l_a_q_fixed_covr <- l_a_q_fixed >= SF[,3] & l_a_q_fixed <= SF[,4]
print(l_a_q_fixed_covr)

## Covr: random effect on l_a_q
SF <- fit$model_fit$summary(variables = 'l_a_q_sigma',
                      "mean",
                      extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
SF
l_a_q_r_sd = simdat$l_a_q_r_sd
l_a_q_r_sd
l_a_q_r_covr <- l_a_q_r_sd >= SF[,3] & l_a_q_r_sd <= SF[,4]
print(l_a_q_r_covr)

## Covr: l_b_q 
SF <- fit$model_fit$summary(variables = 'l_b_q',
                      "mean",
                      extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
SF
l_b_q_fixed = simdat$l_b_q_fixed
l_b_q_fixed
l_b_q_fixed_covr <- l_b_q_fixed >= SF[,3] & l_b_q_fixed <= SF[,4]
print(l_b_q_fixed_covr)

## Covr: random effect on l_a_q
SF <- fit$model_fit$summary(variables = 'l_b_q_sigma',
                      "mean",
                      extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
SF
l_b_q_r_sd = simdat$l_b_q_r_sd
l_b_q_r_sd
l_b_q_r_covr <- l_b_q_r_sd >= SF[,3] & l_b_q_r_sd <= SF[,4]
print(l_b_q_r_covr)

## S
SF <- fit$model_fit$summary(variables = 'S_vec_fixed',
                      "mean",
                      extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
SF
cors <- simdat$Fixed_S[lower.tri(simdat$Fixed_S)]
cors
fix_S <- cors >= SF[,3] & cors <= SF[,4]
print(fix_S )

# Random efx of S
SF <- fit$model_fit$summary(variables = 'S_vec_tau',
                      "mean",
                      extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
SF
ranS_sd = simdat$ranS_sd
ranS_sd
ran_S <- ranS_sd >= SF[,3] & ranS_sd <= SF[,4]
print(ran_S )


###########
## Simulate from fitted stan object
###########

## generated quantities returns num_iteration samples.
## Re-run the the models on the first 100 replications

## Step 1:
## Extract 100 random rows of generated data
rows <- sample(1:1000,  100 )
dt <- data.table(fit$model_fit$draws( variables = "rts_out")[rows,])
dim(dt)

## Position of relevant info:
## [N, tl, nts]

## Eg. first replication 
dt[1, 1]

## init objects:
## list containing all replicaion datasets:
replication_data <- list( )

## objects used to extract the data and assign it to the right unit and variable
tsl <- 1:tl
cols <- 1:nts

## init single replicate data object
result_matrix <- list( )

for(r in 1:100){
  ## Construct a list of length N, containing the tl x nts matrices
  for(i in seq_len(N) ) {
    ## Generate the column names
    column_names <- unlist(lapply(cols, function(col) {
      paste0('rts_out[',i,',', tsl, ',', col, ']')
    }))
    ## Extract the values and transform them into a matrix
    matrix_data <- dt[1, .SD, .SDcols = column_names]
    result_matrix[[i]] <- matrix(unlist(matrix_data), nrow = tl, ncol = nts)
  }
  ## Fill a list that contains all the generated datasets
  replication_data[[r]] = result_matrix 
}

## Sanity check
replication_data[[1]]
replication_data[[2]]



## Init all objects:
phi0_fixed_cov <-phi0_sd_cov <-phi_fixed_covr <-phi_ranef_covr <-log_c_fixed_covr <-log_c_rand_covr <-log_a_fixed_covr <-a_h_r_covr <-log_b_fixed_covr <-b_h_r_covr <-l_a_q_fixed_covr <-l_a_q_r_covr <-l_b_q_fixed_covr <-l_b_q_r_covr <-fix_S <- ran_S <- list( )


## Function to prevent loop to stop when stan model encounters errors
default_return <- NULL
safe_sample <- purrr::possibly( function(s) {
  dcnet( data =  replication_data[[s]], parameterization = 'DCCr' , J =  N,
        group =  groupvec, standardize_data = FALSE, init = 1,
        threads = 1, sampling_algorithm = 'variational')
}, otherwise = default_return)


## Fit model to replication dataset(s)
for(s in 1:100 ) {
  
  fit_r <- safe_sample(s)
  if(!is.null(fit_r$error)) {
    next
  }

 # fit_r = dcnet( data =  replication_data[[s]], parameterization = 'DCCr' , J =  N,
 #               group =  groupvec, standardize_data = FALSE, init = 1,
 #               threads = 1, sampling_algorithm = 'variational'))

              
  ## Coverage: phi0_fixed
  ## Population value(s)
 pop <- fit$model_fit$summary(variables = 'phi0_fixed',
                               "mean",
                               extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  pop
  ## Sampled
  SF <- fit_r$model_fit$summary(variables = 'phi0_fixed',
                                "mean",
                                extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF
  ##
  phi0_fixed <- pop$mean
  phi0_fixed
  phi0_fixed_cov[[s]] <- phi0_fixed >= SF[,3] & phi0_fixed <= SF[,4]
  print(phi0_fixed_cov)
  ##
  ## Coverage: phi0_sd
  pop <- fit$model_fit$summary(variables = 'phi0_tau',
                               "mean",
                               extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  pop
  SF <- fit_r$model_fit$summary(variables = 'phi0_tau',
                                "mean",
                                extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF
  phi0_sd <- pop$mean
  phi0_sd
  phi0_sd_cov[[s]] <- phi0_sd >= SF[,3] & phi0_sd <= SF[,4]
  print(phi0_sd_cov)
  ## Coverage: phi_fixed
  pop <- fit$model_fit$summary(variables = 'vec_phi_fixed',
                               "mean",
                               extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  pop
  SF <- fit_r$model_fit$summary(variables = 'vec_phi_fixed',
                                "mean",
                                extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF
  fixed_phi <- pop$mean
  phi_fixed_covr[[s]] <- fixed_phi >= SF[,3] & fixed_phi <= SF[,4]
  print(phi_fixed_covr)
  ##
  ## Coverage: phi_ranef_sd
  pop <- fit$model_fit$summary(variables = 'phi_tau',
                               "mean",
                               extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF <- fit_r$model_fit$summary(variables = 'phi_tau',
                                "mean",
                                extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF
  phi_ranef_sd = pop$mean
  phi_ranef_sd
  phi_ranef_covr[[s]] <- phi_ranef_sd >= SF[,3] & phi_ranef_sd <= SF[,4]
  print(phi_ranef_covr)
  ##
  ## Covr: Log C fixed: log_c_fixed = rep(0, n_ts), ## on log scale
  pop <- fit$model_fit$summary(variables = 'c_h_fixed',
                               "mean",
                               extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  pop
  SF <- fit_r$model_fit$summary(variables = 'c_h_fixed',
                                "mean",
                                extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF
  log_c_fixed = pop$mean
  log_c_fixed_covr[[s]] <- log_c_fixed >= SF[,3] & log_c_fixed <= SF[,4]
  print(log_c_fixed_covr)
  ##
  ## Covr: random effect on C
  pop <- fit$model_fit$summary(variables = 'c_h_tau',
                               "mean",
                               extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  pop
  SF <- fit_r$model_fit$summary(variables = 'c_h_tau',
                                "mean",
                                extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF
  log_c_r_sd = pop$mean
  log_c_r_sd
  log_c_rand_covr[[s]] <- log_c_r_sd >= SF[,3] & log_c_r_sd <= SF[,4]
  print(log_c_rand_covr)
  ##
  ## Covr: Log A fixed: log_a_fixed = rep(0, n_ts), ## on log scale
  pop <- fit$model_fit$summary(variables = 'a_h_fixed',
                               "mean",
                               extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  pop
  SF <- fit_r$model_fit$summary(variables = 'a_h_fixed',
                                "mean",
                                extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF
  a_h_fixed = pop$mean
  log_a_fixed_covr[[s]] <- a_h_fixed >= SF[,3] & a_h_fixed <= SF[,4]
  print(log_a_fixed_covr)
  ##
  ## Covr: random effect on A
  pop <- fit$model_fit$summary(variables = 'a_h_tau',
                               "mean",
                               extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF <- fit_r$model_fit$summary(variables = 'a_h_tau',
                                "mean",
                                extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF
  a_h_r_sd = pop$mean
  a_h_r_covr[[s]] <- a_h_r_sd >= SF[,3] & a_h_r_sd <= SF[,4]
  print(a_h_r_covr)
  ##
  ## Covr: Log B fixed: log_b_fixed 
  pop <- fit$model_fit$summary(variables = 'b_h_fixed',
                               "mean",
                               extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF <- fit_r$model_fit$summary(variables = 'b_h_fixed',
                                "mean",
                                extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF
  b_h_fixed = pop$mean
  log_b_fixed_covr[[s]] <- b_h_fixed >= SF[,3] & b_h_fixed <= SF[,4]
  print(log_b_fixed_covr)
  ##
  ## Covr: random effect on b
  pop <- fit$model_fit$summary(variables = 'b_h_tau',
                               "mean",
                               extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF <- fit_r$model_fit$summary(variables = 'b_h_tau',
                                "mean",
                                extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF
  b_h_r_sd = pop$mean
  b_h_r_covr[[s]] <- b_h_r_sd >= SF[,3] & b_h_r_sd <= SF[,4]
  print(b_h_r_covr)
  ##
  ## Covr: l_a_q 
  pop <- fit$model_fit$summary(variables = 'l_a_q',
                               "mean",
                               extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF <- fit_r$model_fit$summary(variables = 'l_a_q',
                                "mean",
                                extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF
  l_a_q_fixed = pop$mean
  l_a_q_fixed_covr[[s]] <- l_a_q_fixed >= SF[,3] & l_a_q_fixed <= SF[,4]
  print(l_a_q_fixed_covr)
  ##
  ## Covr: random effect on l_a_q
  pop <- fit$model_fit$summary(variables = 'l_a_q_sigma',
                               "mean",
                               extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF <- fit_r$model_fit$summary(variables = 'l_a_q_sigma',
                                "mean",
                                extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF
  l_a_q_r_sd = pop$mean
  l_a_q_r_covr[[s]] <- l_a_q_r_sd >= SF[,3] & l_a_q_r_sd <= SF[,4]
  print(l_a_q_r_covr)
  ## Covr: l_b_q 
  pop <- fit$model_fit$summary(variables = 'l_b_q',
                               "mean",
                               extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  ##
  SF <- fit_r$model_fit$summary(variables = 'l_b_q',
                                "mean",
                                extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF
  l_b_q_fixed = pop$mean
  l_b_q_fixed_covr[[s]] <- l_b_q_fixed >= SF[,3] & l_b_q_fixed <= SF[,4]
  print(l_b_q_fixed_covr)
  ## Covr: random effect on l_a_q
  pop <- fit$model_fit$summary(variables = 'l_b_q_sigma',
                               "mean",
                               extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF <- fit_r$model_fit$summary(variables = 'l_b_q_sigma',
                                "mean",
                                extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF
  l_b_q_r_sd = pop$mean
  l_b_q_r_covr[[s]] <- l_b_q_r_sd >= SF[,3] & l_b_q_r_sd <= SF[,4]
  print(l_b_q_r_covr)
  ## S
  pop <- fit$model_fit$summary(variables = 'S_vec_fixed',
                               "mean",
                               extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF <- fit_r$model_fit$summary(variables = 'S_vec_fixed',
                                "mean",
                                extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF
  cors <- pop$mean
  cors
  fix_S[[s]] <- cors >= SF[,3] & cors <= SF[,4]
  print(fix_S )
  ## Random efx of S
  pop <- fit$model_fit$summary(variables = 'S_vec_tau',
                               "mean",
                               extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF <- fit_r$model_fit$summary(variables = 'S_vec_tau',
                                "mean",
                                extra_quantiles = ~posterior::quantile2(., probs = c(0.025, 0.975)) )
  SF
  ranS_sd = pop$mean
  ##
  ran_S[[s]] <- ranS_sd >= SF[,3] & ranS_sd <= SF[,4]
  print(ran_S )
}

s


l_a_q_r_covr

## Compute coverage
coverage <- function(param) {
  out <- sapply(1:length(param[[1]]), function(i) {
    mean(sapply(param, function(x) x[i]))
  })
  return(out)
}

covresults <- list()
## 1:
## 2: N = 100, tslength=50
covresults[[2]] <- list(
coverage(phi0_fixed_cov),
coverage(phi0_sd_cov),
coverage(phi_fixed_covr),
coverage(phi_ranef_covr),
coverage(log_c_fixed_covr),
coverage(log_c_rand_covr),
coverage(log_a_fixed_covr),
coverage(a_h_r_covr),
coverage(log_b_fixed_covr),
coverage(b_h_r_covr),
coverage(l_a_q_fixed_covr),
coverage(l_a_q_r_covr),
coverage(l_b_q_fixed_covr),
coverage(l_b_q_r_covr),
coverage(fix_S),
coverage(ran_S))




############################################################
############################################################


fit$model_fit$summary( )



simdat$fixed_phi


## Save out draws and read them back in selectively
fit$model_fit$save_output_files(dir = "../../../Documents/",  basename = "rands",  random = FALSE )

rm("fit")
gc()
## use fread from data.table to handle big files
dt <- fread( input = "../../../Documents/rands-202304061114-1-3bd87d.csv", select = "S_vec_fixed" )

dt <- vroom::vroom( "../../../Documents/rands-202304061114-1-3bd87d.csv")
dt$stan_version_major

dt[1]


varsel <- read_cmdstan_csv( "../../../Documents/rands-202304061114-1-3bd87d.csv" , variables = "S_vec_fixed")


data.table::data.table( fit$model_fit$summary(variables = c('S_vec_fixed',
                                                            'c_h_fixed', 'a_h_fixed', 'b_h_fixed',
                                                            'l_a_q', 'l_b_q')))#, 'Sfixed') ))

w <- ivv(
  tanh(data.table::data.table( fit$model_fit$summary(variables = c('S_vec_fixed') ))$mean)
, nts)
w
w%*%t(w)


lo <- ivv( tanh(data.table::data.table( fit$model_fit$summary(variables = c('S_vec_fixed') ))$q5), nts)
lo%*%t(lo)


hi <- ivv( tanh(data.table::data.table( fit$model_fit$summary(variables = c('S_vec_fixed') ))$q95), nts)
hi%*%t(hi)

vars <- paste0("S_vec_fixed[",  1:6,  "]" )
apply(fit$model_fit$draws( )[, vars], 2, quantile, c(0.025, .975))
hi <- ivv( tanh(apply(fit$model_fit$draws( )[, vars], 2, quantile, c(0.975))), nts)
hi%*%t(hi)


simdat$Fixed_S


quantile(fit$model_fit$draws( )[,'Sfixed[1,2]'], probs = c(.025,  .975))
quantile(fit$model_fit$draws( )[,'Sfixed[1,3]'], probs = c(.025,  .975))
quantile(fit$model_fit$draws( )[,'Sfixed[1,4]'], probs = c(.025,  .975))
quantile(fit$model_fit$draws( )[,'Sfixed[2,3]'], probs = c(.025,  .975))
quantile(fit$model_fit$draws( )[,'Sfixed[2,4]'], probs = c(.025,  .975))
quantile(fit$model_fit$draws( )[,'Sfixed[3,4]'], probs = c(.025,  .975))


quantile(fit$model_fit$draws( )[,'S[1,2]'], probs = c(.025,  .975))
quantile(fit$model_fit$draws( )[,'S[1,3]'], probs = c(.025,  .975))
quantile(fit$model_fit$draws( )[,'S[1,4]'], probs = c(.025,  .975))

colnames(fit$model_fit$draws( ))[grep( "fixed",  colnames(fit$model_fit$draws( )) )]

#Coverage
sum(sign(quantile(fit$model_fit$draws( )[,'Sfixed[1,2]'], probs = c(.025,  .975)) - simdat$Fixed_S[1,2]))== 0
sum(sign(quantile(fit$model_fit$draws( )[,'Sfixed[1,3]'], probs = c(.025,  .975)) - simdat$Fixed_S[1,3]))== 0
sum(sign(quantile(fit$model_fit$draws( )[,'Sfixed[1,4]'], probs = c(.025,  .975)) - simdat$Fixed_S[1,4]))== 0
sum(sign(quantile(fit$model_fit$draws( )[,'Sfixed[2,3]'], probs = c(.025,  .975)) - simdat$Fixed_S[2,3]))== 0
sum(sign(quantile(fit$model_fit$draws( )[,'Sfixed[2,4]'], probs = c(.025,  .975)) - simdat$Fixed_S[2,4]))== 0
sum(sign(quantile(fit$model_fit$draws( )[,'Sfixed[3,4]'], probs = c(.025,  .975)) - simdat$Fixed_S[3,4]))== 0


## Fixed S
sum(sign(quantile(fit$model_fit$draws( )[,'S[1,2]'], probs = c(.025,  .975)) - simdat$Fixed_S[1,2]))== 0
sum(sign(quantile(fit$model_fit$draws( )[,'S[1,3]'], probs = c(.025,  .975)) - simdat$Fixed_S[1,3]))== 0
sum(sign(quantile(fit$model_fit$draws( )[,'S[1,4]'], probs = c(.025,  .975)) - simdat$Fixed_S[1,4]))== 0
sum(sign(quantile(fit$model_fit$draws( )[,'S[2,3]'], probs = c(.025,  .975)) - simdat$Fixed_S[2,3]))== 0
sum(sign(quantile(fit$model_fit$draws( )[,'S[2,4]'], probs = c(.025,  .975)) - simdat$Fixed_S[2,4]))== 0
sum(sign(quantile(fit$model_fit$draws( )[,'S[3,4]'], probs = c(.025,  .975)) - simdat$Fixed_S[3,4]))== 0









## Obtain median R[1,2] across all N
r12 <- NA
for(seqlen in 1:tl ) {
  r12[seqlen] <- median(fit$model_fit$draws( )[, paste0('R[1,', seqlen, ',1,2]')  ])
}
r12

## Compute median pcor across all N from generating simdat
genP <- apply(simdat$DCC_R[1,2,1:tl,],1, median)
genP

genRaw <- apply(simdat$RawCorr[1,2,1:tl,],1, median)
genRaw


## Obtain median R's across all samples and all N from estimated
plout <- sapply(1:N,  function(f) {
  sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('R[',f,',',x,',1,',2,']')]))
})
plout
obsR <- apply(plout, 1, median)

## How close are the estimates?
cor(genP,  obsR, use = "pairwise.complete.obs")

## Transform H to pcor (only for non DIRD
precision12 <- sapply(1:N,  function(f) {
  sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('precision[',f,',',x,',1,',2,']')]))
})
precision11 <- sapply(1:N,  function(f) {
  sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('precision[',f,',',x,',1,',1,']')]))
})
precision22 <- sapply(1:N,  function(f) {
  sapply(seq_len(tl), function(x) median(fit$model_fit$draws( )[, paste0('precision[',f,',',x,',2,',2,']')]))
})

## compute pcor for elements 1,2 
pcor12 <- -apply(precision12, 1, median)/sqrt(apply(precision11, 1, median)*apply(precision22, 1, median))

## How close are the estimates?
cor(genP,  pcor12, use = "pairwise.complete.obs")

cbind(genP, pcor12)
plot(1:tl, genP,  type = 'l', ylim = c(-.3, 0.4))
lines(1:tl, pcor12, col = 'red')


lines(1:tl, obsR, col = 'blue')
lines(1:tl, genRaw, col = 'green')

## Looks like R is indeed tracing the raw correlation among variables and not the pcor.
## pcor12 and genP are indeed pcorrs


median(fit$model_fit$draws( )[, 'S[1,2]'])




cor( cbind(obs = simdat$DCC_R[1,2,1:30,1], gen = r12 )[-1,] )
median(r12)
simdat$S


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

dfs <- read.table("~/UZH/R/Packages/dcnet/tests/local/dataFlip.txt", header = TRUE ) 
head( dfs )

dat1 <- dfs[,c("id","day", "gm_sociable",  "gm_creative", "gm_friendly", "gm_organised")]

## make wide to create NA's and equal length for all
X0 <- reshape(dat1, idvar = "id", timevar = "day", direction = "wide",
        v.names = c("gm_sociable",  "gm_creative", "gm_friendly", "gm_organised"))

varnames = c("gm_sociable",  "gm_creative", "gm_friendly", "gm_organised")
X0
rowSums(is.na(X0 ))
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
max(rowSums( is.na(X2 ) ))
X2
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

head(tsdat[[1]])

plot(tsdat[[50]][,1] , type = 'l')
lines(tsdat[[50]][,2] , col = 'red')

setwd("./tests" )

fit <- dcnet( data = tsdat, J =  N, group =  groupvec, standardize_data = FALSE,
             parameterization = "DCCr", threads = 1, init = 1.5,
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


## Check corr transforms
library(clusterGeneration )
R <- rcorrmatrix(4 )

D <- diag(rgamma(4, 5, 1))
D
H <- D%*%R%*%D
C <- solve(H)
Ds <- diag(sqrt(diag(C)))

-Ds%*%solve(H)%*%Ds+2*diag(rep(1, 4 ))
