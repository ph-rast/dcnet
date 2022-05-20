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

##' @title Simulate TS's for N individuals
##' @param tslength Length of time series, for all N (everyone has same TS length)
##' @param n_ts Number of simultaneous TS per person. Currently only supports 3 or less
##' @param N Subjects
##' @param ranef_sd_S Size of random effects SD in the function that generates random effects around fixed sqrt(r). I'd suggest a rather small number
##' @param phi0_fixed Vector of population values of length n_ts
##' @importFrom clusterGeneration rcorrmatrix

.simuDCC <- function(tslength, n_ts, N, ranef_sd_S, phi0_fixed, alpha =  0.5) {

  ## Define fixed diag for c_h on log scale
  log_c_fixed_diag <- rnorm(n_ts, 0, .5)
  
  ## individual deviations
  log_c_dev_ind <- t(replicate(N, rnorm(n_ts, 0,  .1 ) ))
  
  ## combine to fixed +  ranef
  c_h <- exp( log_c_fixed_diag + log_c_dev_ind )
  
  ## Create individual a_h and b_h
  ab_h_rep <- replicate(N, .ab_h(n_ts ))
  a_h <- ab_h_rep[1,]
  b_h <- ab_h_rep[2,]
  
  ## Same as with a_h and b_h but these ones are for the corrs. Hence, only one value per N
  ab_q_rep <- replicate(N, .ab_h(1))
  a_q <- ab_q_rep[1,]
  b_q <- ab_q_rep[2,]
  
  ## location
  ## Generate random starting values for the location intercept as sum of fixed plus random
  phi0 <- phi0_fixed + replicate(n = N, rnorm(n_ts, 0, .5 ))
  
  ## phi is bound by -1;1. n_tsXn_ts matrix
  ## Create N individual matricesb
  temp <- replicate(n = N, runif(n_ts^2, -.5, .5) )
  phi <- apply(temp, 2, matrix, ncol = n_ts, simplify = FALSE)
  
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
  
  ## Distribution of random effects    
  ## Unconditional Corr;
  ## Fixed :
  Sc <- rcorrmatrix( d =  n_ts, alphad = n_ts)
  
  ## Add random effects with convex method
  
  S <- replicate(N, ( (1-alpha) * Sc + alpha * rcorrmatrix( d =  n_ts, alphad = 100) ) )
  
  ## Old:
  ## Generate N individual unconditional correlation matrices S
  #S <- replicate(N, cov2cor( .covranef(Sc , ranef_sd = ranef_sd_S)) )
    
    ## Q is symmetric
    Q <- array(diag(n_ts), c(n_ts,n_ts, tslength, N))
    Qs <- array(diag(n_ts), c(n_ts,n_ts, tslength, N))
    R <-  array(Sc, c(n_ts,n_ts, tslength, N))

    u <- array(NA, c(n_ts,tslength,N))

    for(j in 1:N ) {
      h[,1,j] <- c_h[j,]
      
      for(t in 2:tslength) {
        DCC_mu[,t,j] <- 
          phi0[,j] + phi[[j]] %*% (y[,t-1,j] - DCC_mu[,t-1,j])

        for(i in 1:n_ts) {
          h[i,t,j] <- sqrt(c_h[j,i] + a_h[[j]][i]*(y[i, t-1,j] - DCC_mu[i, t-1,j])^2 + b_h[[j]][i]*h[i,t-1,j])
        }
        
        u[,t-1,j] <-
          solve( diag(h[,t-1,j]) ) %*% (y[,t-1,j] - DCC_mu[,t-1,j])

        Q[,,t,j] <- (1 - a_q[[j]] - b_q[[j]]) * S[,,j] +
          a_q[[j]] * (u[,t-1,j] %*% t(u[,t-1,j])) + 
          b_q[[j]] * Q[,,t-1,j] 
        
        R[,,t,j] <- cov2cor(Q[,,t,j])

        DCC_H[,,t,j] <- diag(h[,t,j])%*%R[,,t,j] %*%diag(h[,t,j])
        ##
        y[,t,j] <- mvrnorm(mu = DCC_mu[,t,j], Sigma = DCC_H[,,t,j])
        DCC_R[,,t,j] <- cov2cor(DCC_H[,,t,j] )
      }
      DCC_y <- y
    }
    return(list(DCC_y, DCC_R))
  }
