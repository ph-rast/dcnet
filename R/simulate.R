##' Simulate times series
##' @title Create random effects in fixed correlation
##' @param fixed_covmat Fixed covariance matrix, to which ranefs will be added
##' @param ranef_sd SD of random effect that will be added to the cholesky transformed values
##' @author Philippe Rast
##' @importFrom ks vech 
##' @importFrom ks invvech


.covranef <- function(fixed_covmat,  ranef_sd ) {
    c_f <- fixed_covmat
    ## Add individual deviations from c_f:
    ## Its probably best to cholesky factor, add ranef deviations, an recompute cov:
    
    ## vectorize cholesky of c_f, add random effects with rnorm, save as Lv
    Lv <- vech(t(chol(c_f ))) + rnorm( length(vech(t(chol(c_f )))), 0, ranef_sd )
    ## Recompute cholesy as matrix
    Lvi <- invvech( Lv )
    ## invvech also returns values in upper triangle; set to zero
    Lvi[upper.tri(c_f )] <- 0
    ## Compute covariance matrix
    sigma <- Lvi %*% t(Lvi)
    return( sigma )
}

##' Draw global and local scales for the (half‑)Student‑t / Cauchy block horseshoe
##' @param p probability for rcauchy and rt
##' @param global_df global df's
##' @param local_df local df
##' @param g_scale global scale
##' @param l_scale local scale
.hs_scales <- function(p,
                       global_df = 1,   # df = 1 for half‑Cauchy
                       local_df  = 1,
                       g_scale   = 0.1, # scale of the global tau
                       l_scale   = 1) { # scale of local lambdas
    ## global scale tau
    tau <- if (global_df == 1) {
               abs(rcauchy(1, location = 0, scale = g_scale))
           } else {
               abs(rt(1, df = global_df) * g_scale)
           }
    ## local scales gammas
    lambda <- if (local_df == 1) {
                  abs(rcauchy(p, location = 0, scale = l_scale))
              } else {
                  abs(rt(p, df = local_df) * l_scale)
              }

    list(tau = tau, lambda = lambda)
}

##' Create Random correlations
##' @param n_ts Number of time series
##' @author Philippe Rast
.ranS <- function(n_ts, fixedS, ranS_sd) {
    lt <- choose(n_ts, 2)
    x <- tanh(atanh(fixedS) + rnorm(lt, 0, sd = ranS_sd))
    z <- diag(rep(0, n_ts))
    index <- 0
    for (j in 1:n_ts) {
        for (i in 1:n_ts) {
            if (i > j) {
                index <- index + 1
                z[i, j] <- x[index]
            }
            if (i < j) {
                z[i, j] <- 0
            }
        }
    }
    L <- diag(rep(0, n_ts))
    for (i in 1:n_ts) {
        for (j in 1:n_ts) {
            if (i < j) {
                L[i, j] <- 0.0
            }
            if (i == j) {
                if (i == 1) {
                    L[i, j] <- 1.0
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
    R <- L %*% t(L)
    return(R)
}

##' @title Internal function to generate constrained a's and b's
##' @description Function to generate individual a_h and b_h of lenght n_ts that comply with a_h + b_h < 1 
##' @param n_ts Number of time series
##' @importFrom MASS mvrnorm
##' @author Philippe Rast
.ab_h <- function(n_ts ) {
    a_h <- runif(n_ts, 0, 1 )
    ## Ensure that sum of a_h and b_h is !> 1
    b_h <- runif(n_ts, 0, 1-a_h )
    return(list(a_h, b_h ) )
}

##' @title Create phi matrix that ensures stationarity of VAR(1)
##' The function returns a matrix with eigen values less than 1 in modulus
##' @param fixed_phi A fixed effects square matrix
.stat_var <- function(n_ts, fixed_phi, phi_ranef_sd, D, stationarity = TRUE) {
    ##D <- diag(tanh(rnorm(n_ts, 0, 1) ) )
    vecfix <- c(fixed_phi) ## columnwise
    indfix <- vecfix + rnorm(n_ts^2, 0, phi_ranef_sd)
    A <- matrix(indfix, nrow = n_ts)
    stationary <- A %*% D %*% solve(A)
    nonstationary <- A
    if (stationarity == TRUE) {
        return(stationary)
    } else return(nonstationary)
}


##' @title Simulate TS's for N individuals
##' @param tslength Length of time series, for all N (everyone has same TS length)
##' @param n_ts Number of simultaneous TS per person. Currently only supports 3 or less
##' @param N Subjects
##' @param ranS_sd Size of random effects SD in the function that generates random effects around fixed sqrt(r). I'd suggest a rather small number
##' @param phi0_fixed Vector of population values of length n_ts
##' @param empirical Sample from the emprirical multivariate normal distribution. Defaults to `FALSE`
##' @importFrom clusterGeneration rcorrmatrix

.simuDCC <- function(tslength, n_ts, N,
                     phi0_fixed = rep(0, n_ts),
                     phi0_sd = .5,
                     phi_mu = 0, ## Values of phi
                     phi_sd_diag =  0.3, ## fixed-effects spread of the diagonal own-lag
                     phi_sd_off = 0.05, ## fixed-effects spread of the cros-lag
                     phi_ranef_sd = 0.01, ## This is the random effect of phi
                     stationarity_phi = FALSE, ## Should the phi values ensure that ts is stationary
                     log_c_fixed = rep(0, n_ts), ## on log scale
                     log_c_r_sd = 0.1,
                     a_h_fixed = rep(0, n_ts),
                     a_h_r_sd = 0.1, ## on logit scale
                     b_h_fixed = rep(0, n_ts),
                     b_h_r_sd = 0.1, ## on logit scale
                     l_a_q_fixed = 0,
                     l_a_q_r_sd = 0.1,
                     l_b_q_fixed = 0,
                     l_b_q_r_sd = 0.1,
                     ranS_sd = 1,
                     empirical = FALSE) {
    ## Define fixed diag for c_h on log scale
    log_c_fixed_diag <- log_c_fixed
    ## individual deviations
    log_c_dev_ind <- t(replicate(N, rnorm(n_ts, 0,  log_c_r_sd)))
    ## combine to fixed +  ranef
    c_h <- exp(log_c_fixed_diag + log_c_dev_ind)

    ## Create individual a_h and b_h
    a_h_random <- t(replicate(N, rnorm(n_ts, 0,  a_h_r_sd)))
    a_h <- 1 / (1 + exp(- (a_h_fixed + a_h_random)))

    ## Create individual a_h and b_h
    b_h_random <- t(replicate(N, rnorm(n_ts, 0,  b_h_r_sd)))
    b_h <- (1 - a_h) / (1 + exp(- (b_h_fixed + b_h_random)))

    a_q <- 1 / (1 + exp(- (l_a_q_fixed + rnorm(N, 0, l_a_q_r_sd))))
    b_q <- (1 - a_q) / (1 + exp(- (l_b_q_fixed + rnorm(N, 0, l_b_q_r_sd))))

    ## location
    ## Generate random starting values for the location intercept as sum of fixed plus random
    phi0 <- phi0_fixed + replicate(n = N, rnorm(n_ts, 0, phi0_sd))

    ## phi is bound by -1;1. n_tsXn_ts matrix
    ## Create N individual matrices
    ## watch out, this might create non-stationary series
    ## To avoid this, create matrices who's eigenvalues are less than 1 in modulus.
    ## This ensures stationary of the VAR(1) model: https://phdinds-aim.github.io/time_series_handbook/03_VectorAutoregressiveModels/03_VectorAutoregressiveMethods.html#stationarity-of-the-var-1-model

    ## Create fixed phi
    ## phi_fixed <- matrix(rnorm(n_ts^2, phi_mu, phi_sd), nrow = n_ts)
    ## --- NEW: diagonal‑vs‑off‑diagonal prior ------------------------------------
    ## user‑tunable hyper‑SDs (defaults give old behaviour w.out off-diagonal shrinkage)
    phi_sd_diag <- phi_sd_diag # spread of own‑lag coefficients
    phi_sd_off <- phi_sd_off # baseline spread of cross‑lags *before* horseshoe

    ## sample diagonal (own‑lags) – usually positive and comparatively large
    diag_vals <- rnorm(n_ts, mean = phi_mu, sd = phi_sd_diag)

    ## sample block‑horseshoe scales for n_ts*(n_ts‑1) off‑diagonals
    hs <- .hs_scales(p = n_ts^2 - n_ts, g_scale = phi_sd_off)
    off_vals <- rnorm(n_ts^2 - n_ts, 0, sd = hs$tau * hs$lambda)

    ## assemble matrix
    phi_fixed <- matrix(0, n_ts, n_ts)
    diag(phi_fixed) <- diag_vals
    phi_fixed[upper.tri(phi_fixed) | lower.tri(phi_fixed)] <- off_vals
    ## ---------------------------------------------------------------------------
    
    D <- diag(tanh(rnorm(n_ts, 0, 1)))

    ## add  random effect with .stat_var function
    phi <- replicate(N, .stat_var(n_ts,
                                  fixed_phi = phi_fixed,
                                  phi_ranef_sd = phi_ranef_sd,
                                  D,
                                  stationarity = stationarity_phi),
                     simplify = FALSE)

  y <- array(NA, dim = c(n_ts, tslength, N))
  ## create named variables:
  rownames(y) <-  paste0('X',  seq_len(n_ts) )

  ## init y
  y[,1,] <- phi0
  
  ## mean
  DCC_mu <- array(0, dim = c(n_ts, tslength, N))
  
  h <- array(.5, dim = c(n_ts, tslength, N))
  DCC_H <- array( NA, c(n_ts,n_ts, tslength, N))
  DCC_R <- array( NA, c(n_ts,n_ts, tslength, N))
  RawCorr <- array( NA, c(n_ts,n_ts, tslength, N))
  
    ## Distribution of random effects
    ## Unconditional Corr;
    ## Fixed :
    Sc <- clusterGeneration::rcorrmatrix( d =  n_ts, alphad = n_ts*2)
    Fixed_S <- Sc[lower.tri(Sc)]

    S <- list()
    for (j in 1:N) {
        S[[j]] <- .ranS(n_ts = n_ts, fixedS = Fixed_S, ranS_sd = ranS_sd)
    }

    ## Q is symmetric
    Q <- array(diag(n_ts), c(n_ts, n_ts, tslength, N))
    Qs <- array(diag(n_ts), c(n_ts, n_ts, tslength, N))
    R <-  array(S[[1]], c(n_ts, n_ts, tslength, N))

    u <- array(NA, c(n_ts, tslength,N))

    for (j in 1:N) {
        h[, 1, j] <- c_h[j, ]
      
        for (t in 2:tslength) {
            DCC_mu[, t, j] <-
                phi0[, j] + phi[[j]] %*% (y[, t - 1, j] - phi0[, j]) #DCC_mu[, t - 1, j])

            for (i in 1:n_ts) {
                h[i,t,j] <-
                    sqrt(c_h[j,i] +
                         a_h[j,i] * (y[i, t - 1, j] - DCC_mu[i, t - 1, j])^2 +
                         b_h[j,i] * h[i, t - 1, j]^2)
            }

            u[, t - 1, j] <-
                solve(diag(h[, t - 1, j])) %*% (y[, t - 1, j] - DCC_mu[, t - 1, j])

            Q[, , t, j] <-
                (1 - a_q[[j]] - b_q[[j]]) * S[[j]] +
                a_q[[j]] * (u[, t - 1, j] %*% t(u[, t - 1, j])) +
              b_q[[j]] * Q[, , t - 1, j]

            R[, , t, j] <- cov2cor(Q[, , t, j])
            R[, , t, j]

            DCC_H[, , t, j] <-
                diag(h[, t, j]) %*% R[, , t, j] %*% diag(h[, t, j])
            DCC_H[, , t, j]
            ##
        y[,t,j] <- mvrnorm(mu = DCC_mu[,t,j], Sigma = DCC_H[,,t,j], empirical =  empirical)
        ## compute partial corr from DCC_H
        DCC_R[,,t,j] <- cov2cor( solve(DCC_H[,,t,j]) ) * ( diag(n_ts) - 1)
        ## Obtain simple corr
        RawCorr[,,t,j] <- cov2cor(DCC_H[,,t,j])
      }
      DCC_y <- y
    }
  if(stationarity_phi == TRUE ) {
    fixed_phi <- phi_fixed # %*%D%*%solve(phi_fixed)
  } else {
    fixed_phi <- phi_fixed
  }
    ## Return lower tri of Fixed_S on tanh scale, for comparison w. model
    fixed_S_atanh <- atanh(Sc[lower.tri(Sc)])
    
  return(list(DCC_y,
              DCC_R = DCC_R,
              S = S,
              RawCorr =  RawCorr,
              Fixed_S = Sc,
              fixed_S_atanh = fixed_S_atanh,
              fixed_phi = fixed_phi,
              ranS_sd = ranS_sd,
              phi_sd_diag =  phi_sd_diag,
              phi_sd_off =  phi_sd_off,
              l_b_q_r_sd = l_b_q_r_sd,
              l_b_q_fixed = l_b_q_fixed,
              l_a_q_r_sd = l_a_q_r_sd,
              l_a_q_fixed = l_a_q_fixed,
              b_h_r_sd = b_h_r_sd,
              b_h_fixed = b_h_fixed,
              a_h_r_sd = a_h_r_sd,
              a_h_fixed = a_h_fixed,
              log_c_r_sd = log_c_r_sd,
              log_c_fixed = log_c_fixed,
              phi_ranef_sd = phi_ranef_sd,
              phi0_sd = phi0_sd,
              phi0_fixed = phi0_fixed))
}
