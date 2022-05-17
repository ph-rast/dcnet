#' Create dcnet data object. Data is standardized (default) across individual variables to facilitate computation.
#' Data object results in a list. 
#' 
#' @param data Multilevel time-series data in long format
#' @param J numeric. Number of groups (individuals)
#' @param group Vector with group id of full length
#' @param xC Numeric matrix containing predictors.
#' @param P Numeric. Currently fixed to 1
#' @param Q Numeric. Currently fixed to 1
#' @param standardize_data Logical.
#' @param distribution Numeric: 0 = Gaussian, 1 = Student-T
#' @param meanstructure Character. Default is multilevel 'VAR'.
#' @param simplify_ch Numeric value indicating whether random effects correlations should be simplfied to just diagonal matrix: 1 = yes; 0 = no.
#' @param simplify_ah Numeric value indicating whether random effects correlations should be simplfied to just diagonal matrix: 1 = yes; 0 = no.
#' @param simplify_bh Numeric value indicating whether random effects correlations should be simplfied to just diagonal matrix: 1 = yes; 0 = no.
#' @return dcnet stan data list. 
#' @keywords internal
stan_data <- function(data, J, group, xC, P = 1, Q = 1, standardize_data, distribution = 0,
                      meanstructure =  'VAR', simplify_ch = 1, simplify_ah = 1, simplify_bh = 1) {

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
    xC <- matrix(0, ncol = nt, nrow = J * tl)
  } else if ( !is.null(xC) ) {
      ## Test dimension of XC
      if( ncol(xC) != nt )  stop("xC must have same dimension as data (same number of columns)")
      if( nrow(xC) != J*tl) stop("xC must have same dimension as data (same number of rows)")
    }
  ## Match dimension of predictor to TS. If only one vector is given, it's assumed that it is the same for all TS's
  ## if ( is.null(ncol(xC)) ) {
  ##   warning("xC is assumed constant across TS's")
  ##   xC <- matrix(xC, nrow = nrow(data), ncol = ncol(data)) ## Risky, better to throw an error
  ## } else if ( dim( xC )[2] != dim( data )[2] ) { ## xC is not a vector  - check if it is of right dimension
  ##   warning("xC is not of right dimension - adapt xC dimension to match number of TS")
  ## }

  ## TODO: Generalize to unequal numbers of observations per individual
  ## WARNIGN: This assumes equal numbers of observations for everyone!
  T <- nrow(data[[1]])
  nobs <- T * J

  ## Check whether all J have equal nrow
  varT <- NA
  for(i in 1:J ) {
    varT[i] <- nrow(data[[i]])
  }
  if( var(varT) > 0) stop("Time-series lengths are not equal for all groups/subjects")
  
  if( standardize_data ) {
    ## Standardize time-series
    stdx <- data
    centered_data <- matrix(NA, nrow = J, ncol = nt)
    scaled_data <- matrix(NA, nrow = J, ncol = nt)
    
    for( i in 1:J ) {
      stdx[[i]] <- scale(data[[i]] )
      centered_data[i,] <- attr(stdx[[i]], "scaled:center")
      scaled_data[i,] <- attr(stdx[[i]], "scaled:scale")  
    }
    
    return_standat <- list(T = T, ## TODO: may need to allow different T's per J
                           rts = stdx,
                           xC = xC,
                           nt = ncol(stdx[[1]]),
                           centered_data = centered_data,
                           scaled_data = scaled_data,
                           distribution = distribution,
                           P = P,
                           Q = Q,
                           meanstructure = meanstructure,
                           J = J,
                           nobs = nobs,
                           group = group,
                           simplify_ch = simplify_ch,
                           simplify_ah = simplify_ah,
                           simplify_bh = simplify_bh)
  } else {
    ## Unstandardized
    return_standat <- list(T = T, ## TODO: may need to allow different T's per J
                           rts = data,
                           xC = xC,
                           nt = ncol(data[[1]]),
                           distribution = distribution,
                           P = P,
                           Q = Q,
                           meanstructure = meanstructure,
                           J = J,
                           nobs = nobs,
                           group = group,
                           simplify_ch = simplify_ch,
                           simplify_ah = simplify_ah,
                           simplify_bh = simplify_bh)
  }
  
  return(return_standat)
}