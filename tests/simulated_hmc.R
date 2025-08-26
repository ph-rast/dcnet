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
N <- 50
tl <- 50
nts <- 3
simdat <- .simuDCC(
    tslength = tl, N = N, n_ts = nts,
    phi_mu = 0, ## populate phi
    phi0_sd = 0.4, ## create variation in the intercepts of the time series
    phi_sd_diag = 0.3, ## SD of own lags across TS's
    phi_sd_off = 0.05, ## SD of cross lags across TS's
    phi_ranef_sd = 0.01, ## random effects in phi
    log_c_fixed = rep(-0.9, nts),
    log_c_r_sd = 0.3,
    a_h_fixed = rep(-2.5, nts),
    a_h_r_sd = 0.1,
    b_h_fixed = rep(-1.5, nts), ## On logit scale
    b_h_r_sd = 0.1,
    l_a_q_fixed = -1.5, ## on logit scale
    l_b_q_fixed = -.5, ## on logit scale
    l_a_q_r_sd = 0.2,
    l_b_q_r_sd = 0.2,
    phi0_fixed = rep(0, nts),
    ranS_sd = 0.25, ## random effects on atanh scale
    stationarity_phi = FALSE)

rtsgen <- lapply(seq(dim(simdat[[1]])[3]), function(x) t(simdat[[1]][, , x]))

groupvec <- rep(c(1:N),  each = tl)


range(lapply(1:N, function(x) range(rtsgen[[x]])))


fit0 <- dcnet(
    data = rtsgen, parameterization = "DCCms", J =  N,
    group = groupvec,
    init = 0.5,
    meanstructure = "VAR",
    iterations = 30000,
    chains = 4,
    threads = 1, # tol_rel_obj =  0.01, ## 8 threads: 188 mins /
    sampling_algorithm = "variational", # "pathfinder", # "laplace",
    grainsize = 1)



summary(fit0)

#### Garch h params
draws_df <- fit0$model_fit$draws(format = "draws_df")


######################################
## Simulate 
######################################

## list containing all replicaion datasets:
replication_data <- list()

## Create data:
simulate_data <- function(N = 100, tl = 50, nts = 3) {
    simdat <- .simuDCC(
        tslength = tl, N = N, n_ts = nts,
        phi_mu = 0, ## populate phi
        phi0_sd = 0.4, ## create variation in the intercepts of the time series
        phi_sd_diag = 0.3,  ## SD of own lags across TS's
        phi_sd_off = 0.05,  ## SD of cross lags across TS's
        phi_ranef_sd = 0.01, ## random effects in phi
        log_c_fixed = rep(-0.9, nts),
        log_c_r_sd = 0.3,
        a_h_fixed = rep(-2.5, nts),
        a_h_r_sd = 0.1,
        b_h_fixed = rep(-1.5, nts),  ## On logit scale
        b_h_r_sd = 0.1,
        l_a_q_fixed = -1.5,  ## on logit scale
        l_b_q_fixed = -.5,   ## on logit scale
        l_a_q_r_sd = 0.2,
        l_b_q_r_sd = 0.2,
        phi0_fixed =  rep(0, nts),
        ranS_sd = 0.25, ## random effects on atanh scale
        stationarity_phi = FALSE
    )
    rtsgen <- lapply(seq(dim(simdat[[1]])[3]), function(x) t(simdat[[1]][, , x]))
    groupvec <- rep(c(1:N), each = tl)
    return(list(rtsgen, groupvec, simdat, N = N))
}


replication_data <- simulate_data()
replication_data[[1]]

replication_data[[2]]
replication_data[[3]]
replication_data$N

## Function to prevent loop to stop when stan model encounters errors
## Try fitting model 3 times before moving on

## V2:
safe_sample <- function(s, max_retries = 1, replication_data) {
    replication_data <- replication_data

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
            dcnet(
                data = replication_data[[1]], parameterization = "DCCms", J = replication_data$N,
                group = replication_data[[2]], standardize_data = FALSE,
                init = 0,
                meanstructure = "VAR",
                iterations = 500,
                #eta = 0.1,
                #tol_rel_obj = 0.007,
                sampling_algorithm = "hmc"
            )
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
      #      "All proposed step-sizes failed",
      #      "algorithm may not have converged",
      #      "Exception:.*not positive definite",
      #      "Exception:",
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
## phi0_fixed_cov <- phi0_sd_cov <-phi_fixed_covr <-phi_ranefdiag_covr<-phi_ranefoff_covr <-log_c_fixed_covr <-log_c_rand_covr <-log_a_fixed_covr <-a_h_r_covr <-log_b_fixed_covr <-b_h_r_covr <-l_a_q_fixed_covr <-l_a_q_r_covr <-l_b_q_fixed_covr <-l_b_q_r_covr <-fix_S <- ran_S <- list()



## Function to compute coverage probability ranging from L to U in the original population distribution
overlap <- function(population, L, U) {
    coverage <- as.numeric(population > L & population < U)
    return(coverage)
}

crirange <- function(L, U) {
    bwidth <- abs(as.numeric(U-L))
    return(bwidth)
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

looic <- function(fit) {
    log_lik_r <- fit$model_fit$draws(variable = "log_lik")
    r_eff_r <- loo::relative_eff(exp(log_lik_r ))
    fr <- loo::loo(log_lik_r, r_eff = r_eff_r)
    looic <- fr$estimates[3,1]
    return(looic)
}

## Stan variables:
variables_m <- c(
    'phi0_fixed', 'phi0_tau', 'vec_phi_fixed', 'sigma_re_own', 'sigma_re_cross',
    'tau_own', 'tau_cross' )
   # 'c_h_fixed', 'c_h_tau', 'a_h_fixed', 'a_h_tau', 'b_h_fixed', 'b_h_tau',
   # 'l_a_q', 'l_a_q_sigma', 'l_b_q', 'l_b_q_sigma', 'S_vec_fixed', 'S_vec_tau')

## Simulation variables: Note that the simulatin script only defines one random
## effect for phi, phi_ranef_sd, but the stan model captures the random effects
## in both, the own and cross lag of phi in sigma_re_own/cross
var <- c(
    "phi0_fixed", "phi0_sd", "fixed_phi", "phi_ranef_sd", "phi_ranef_sd",
    "phi_sd_diag", "phi_sd_off" )
    #"log_c_fixed", "log_c_r_sd", "a_h_fixed", "a_h_r_sd", "b_h_fixed", "b_h_r_sd",
    #"l_a_q_fixed", "l_a_q_r_sd", "l_b_q_fixed", "l_b_q_r_sd", "fixed_S_atanh", "ranS_sd")

if(!length(variables_m) == length(var)) stop("variable list does not match!")

# Initialize empty lists to store the results for each variable
cov_list <- list()
crir_list <- list()
rmse_list <- list()
bias_list <- list()
bins_list <- list()
looic_list <- list()


for (s in 1:10) {

    replication_data <- simulate_data(N = 50, tl = 50)
    fit_r <- safe_sample(s, replication_data = replication_data)
    
    if (is.null(fit_r)) {
        next
    }
    
    for (p in 1:length(variables_m)) {
        ## Add logic for S_vec_tau: Needs to be taken from step 1 fit, all others form step 2
        if (variables_m[p] == "S_vec_tau") {
            SF <- data.frame(t(apply(
                fit_r$S_vec_tau_post, 2,
                function(x) quantile(x, c(0.5, 0.025, .975))
            )))
            colnames(SF) <- c("mean", "q2.5", "q97.5")
        } else {
            SF <- fit_r$model_fit$summary(
                                      variables = variables_m[p], "mean",
                                      extra_quantiles = ~ posterior::quantile2(., probs = c(0.025, 0.975))
                                  )
        }
        cov_list[[variables_m[p]]][[s]] <-
            overlap(
                unlist(replication_data[[3]][var[p]]),
                L = SF$q2.5, U = SF$q97.5
            )
        crir_list[[variables_m[p]]][[s]] <-
            crirange(L = SF$q2.5, U = SF$q97.5)
        rmse_list[[variables_m[p]]][[s]] <-
            sapply(seq_len(nrow(SF)), function(i) {
                rmse(
                    fit_r$model_fit$draws(variables = variables_m[p])[, , i],
                    unlist(replication_data[[3]][var[p]])[i]
                )
            })
        bias_list[[variables_m[p]]][[s]] <-
            sapply(seq_len(nrow(SF)), function(i) {
                bias(
                    fit_r$model_fit$draws(variables = variables_m[p])[, , i],
                    unlist(replication_data[[3]][var[p]])[i]
                )
            })
        looic_list[[variables_m[p]]][[s]] <-
            suppressWarnings(looic(fit_r))
    }
    gc()
    s
}

## Try reducing the convergence criterion to 0.005
cov_list
crir_list
rmse_list
bias_list
bins_list
looic_list

var_averages <- list()
for (i in 1:length(cov_list)) {
    nested_average <- c()                                             
    for (j in 1:length(cov_list[[i]])) {
        if (is.null(cov_list[[i]][[j]])) cov_list[[i]][[j]] <- NA
        nested_average <- c(nested_average, mean(cov_list[[i]][[j]], na.rm = TRUE))
        var_averages[[i]] <- nested_average
    }
}
names(var_averages) <- names(cov_list)
var_averages


var_av <- sapply(var_averages, function(x ) mean(x, na.rm = TRUE))
var_av


bias_averages <- list()
for (i in 1:length(bias_list)) {
    nested_average <- c()
    for (j in 1:length(bias_list[[i]])) {
        if (is.null(bias_list[[i]][[j]])) bias_list[[i]][[j]] <- NA
        nested_average <- c(nested_average, mean(bias_list[[i]][[j]], na.rm = TRUE))
        bias_averages[[i]] <- nested_average
    }
}
names(bias_averages) <- names(bias_list)
bias_averages

bias_av <- sapply(bias_averages, function(x) mean(x, na.rm = TRUE))
round(cbind(bias_av, var_av), 3)

n30tl50 <- round(cbind(bias_av, var_av), 3)
n30tl50
##n50tl50 <- round(var_av, 2)
n50tl50
## n50tl75 <- round(var_av, 2)
n50tl75
## n75tl75 <- round(var_av, 2)
n75tl75
## n75tl100 <- round(var_av, 2)
n75tl100
n100tl100 <- round(cbind(bias_av, var_av), 3)
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


