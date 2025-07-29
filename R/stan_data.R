#' Create dcnet data object. Data is standardized (default) across individual variables to facilitate computation.
#' Data object results in a list. 
#' 
#' @param data Multilevel time-series data in long format
#' @param J numeric. Number of groups (individuals)
#' @param group Vector with group id of full length
#' @param xC Numeric matrix containing predictors.
#' @param S_pred Integer matrix containing predictors for each datapoing.
#' @param P Numeric. Currently fixed to 1
#' @param Q Numeric. Currently fixed to 1
#' @param standardize_data Logical.
#' @param distribution Numeric: 0 = Gaussian, 1 = Student-T
#' @param meanstructure Character. Default is multilevel 'VAR'.
#' @param simplify_ch Numeric value indicating whether random effects correlations should be simplfied to just diagonal matrix: 1 = yes; 0 = no.
#' @param simplify_ah Numeric value indicating whether random effects correlations should be simplfied to just diagonal matrix: 1 = yes; 0 = no.
#' @param simplify_bh Numeric value indicating whether random effects correlations should be simplfied to just diagonal matrix: 1 = yes; 0 = no.
#' @param lbound Lower bound on observed values (across data matrix) - defaults to min(data).
#' @param ubound Upper bound on observed values (across data matrix) - defaults to max(data).
#' @return dcnet stan data list. 
#' @keywords internal
stan_data <- function(data, J, group, xC, S_pred, P = 1, Q = 1, standardize_data, distribution = 0,
                      meanstructure =  'VAR', simplify_ch = 1, simplify_ah = 1, simplify_bh = 1,
                      lbound, ubound, grainsize) {

  ## Check for data type:
  ## dim data needs to return NULL as the list is a collection of individual matrices for each individual.
  if( typeof(data) != "list" | !is.null( dim(data) ) ) {
    stop("Data object is not a list method object.")
  }

  ## Check if content is numeric:
  
  ## Check whether var(variables)>0:
  
  ##  if ( is.null( colnames( data ) ) ) colnames( data ) = paste0('t', 1:ncol( data ) )

  ## Model for meanstructure
  if( meanstructure == "constant" | meanstructure == 0 ) {
    meanstructure <- 0
  } else if ( meanstructure == "arma" | meanstructure == "ARMA" | meanstructure == 1 ){
    meanstructure <- 1
  } else if ( meanstructure == "var"  | meanstructure == "VAR"  | meanstructure == 2) {
    meanstructure <- 2
  } else {
    stop("meanstructure must be either 'constant', 'ARMA' or 'VAR'.")
  }
  
  
  ## Tests on predictor
  nt <- ncol(data[[1]])
    
  ## Pass in a 0 matrix, so that stan does not complain
  
  if ( is.null(xC) ) {
    xC <- matrix(0, ncol = nt, nrow = J * nrow(data[[1]]))
  } else if ( !is.null(xC) ) {
      ## Test dimension of XC
      if( ncol(xC) != nt )  stop("xC must have same dimension as data (same number of columns)")
    if( nrow(xC) != J*nrow(data[[1]])) stop("xC must have same dimension as data (same number of rows)")
  }

  if ( is.null(S_pred) ) {
    S_pred <- lapply(1:J, FUN =  function(x) rep(0, nrow(data[[1]])))
  }
  
  ## Match dimension of predictor to TS. If only one vector is given, it's assumed that it is the same for all TS's
  ## if ( is.null(ncol(xC)) ) {
  ##   warning("xC is assumed constant across TS's")
  ##   xC <- matrix(xC, nrow = nrow(data), ncol = ncol(data)) ## Risky, better to throw an error
  ## } else if ( dim( xC )[2] != dim( data )[2] ) { ## xC is not a vector  - check if it is of right dimension
  ##   warning("xC is not of right dimension - adapt xC dimension to match number of TS")
  ## }

  ## TODO: Generalize to unequal numbers of observations per individual
  ## WARNING: This assumes equal numbers of observations for everyone!
  T <- nrow(data[[1]])
  nobs <- T * J

  ## Check whether all J have equal nrow
  varT <- NA
  for(i in 1:J ) {
    varT[i] <- nrow(data[[i]])
  }
  if( var(varT) > 0) stop("Time-series lengths are not equal for all groups/subjects")

  ## Lower and upper bound on observed data matrix: Min and max across all columns
  if( !lbound ) {
    lbound <- min( unlist(data) )
  }
  if( !ubound ) {
    ubound <- max( unlist(data) )
  }  
  ## Standardization and Centering: Maintain relative locations of time series across individuals/groups to
  ## obtain meaningful random intercept terms.
  ## Same for centering only (when not standardizing). Relative position needs to be maintained.
  ## Obtain Means and SD across J (grouping) within the same variable. 
  if( standardize_data ) {
    ## Standardize time-series. 
    stdx <- data

    ##  Obtain column means across all J
    col_means <- colMeans( do.call(rbind,  stdx) )
    sd_means <- apply( do.call(rbind,  stdx), 2, sd )

    ## First convert integers to numeric, then standardize    
    for( i in 1:J ) {
      stdx[[i]] <- apply( data[[i]], MARGIN = 2, as.numeric ) ## make numeric
      stdx[[i]] <- ( stdx[[i]] - col_means[col(stdx[[i]])] ) / sd_means[col(stdx[[i]])]  ## Vector col_means needs to be subtracted row-wise and divided row_wise
    }

    ## Scale bounds
    lbound <- min( unlist(stdx) )
    ubound <- max( unlist(stdx) )
    
    return_standat <- list(T = T, ## TODO: may need to allow different T's per J
                           rts = stdx,
                           xC = xC,
                           S_pred = S_pred,
                           nt = ncol(stdx[[1]]),
                           grand_mean = col_means, 
                           grand_sd = sd_means, 
                           distribution = distribution,
                           P = P,
                           Q = Q,
                           meanstructure = meanstructure,
                           J = J,
                           nobs = nobs,
                           group = group,
                           simplify_ch = simplify_ch,
                           simplify_ah = simplify_ah,
                           simplify_bh = simplify_bh,
                           lbound = lbound, ubound = ubound,
                           grainsize = grainsize)
  } else {
    ## Unstandardized
    ctrx <- data

    ##  Obtain column means across all J
    col_means <- colMeans( do.call(rbind,  ctrx) )
    
    ## First convert integers to numeric, then subtract grand col_means to center
     for( i in 1:J ) {
       ctrx[[i]] <- apply( data[[i]], MARGIN = 2, as.numeric ) #scale(data[[i]], center = FALSE, scale = FALSE)
       ## No need to center: Stan uses first value of TS as prior for mean
       ## ctrx[[i]] <- ctrx[[i]] - col_means[col(ctrx[[i]])]  ## Vector col_means needs to be subtracted row-wise 
     }
    
    return_standat <- list(T = T, ## TODO: may need to allow different T's per J
                           rts = ctrx,
                           xC = xC,
                           S_pred = S_pred,
                           nt = ncol(ctrx[[1]]),
                           grand_mean = col_means,
                           distribution = distribution,
                           P = P,
                           Q = Q,
                           meanstructure = meanstructure,
                           J = J,
                           nobs = nobs,
                           group = group,
                           simplify_ch = simplify_ch,
                           simplify_ah = simplify_ah,
                           simplify_bh = simplify_bh,
                           lbound = lbound, ubound = ubound,
                           grainsize = grainsize)
  }
  
  return(return_standat)
}
