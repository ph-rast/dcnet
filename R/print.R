##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param object 
##' @return 
##' @author Philippe Rast
##' @import data.table
summary.dcnet <- function(object, CrI = c(.025, .975), digits = 2,  ... ) {
  
  ## Parameters for each model
  common_params <- c("lp")
  var_params <- c("phi0_fixed|phi0_L|phi0_tau|vec_phi_fixed|phi_L|phi_tau")

  ccc_params <- c("c_h|a_h|b_h|R|beta|c_h_var")
  dcc_D_params <- c("c_h_fixed|a_h_fixed|b_h_fixed|c_h_L|c_h_tau|a_h_L|a_h_tau|b_h_L|b_h_tau")
  dcc_Q_params <- c("l_a_q|l_a_q_r|l_b_q|l_b_q_r|S")
  dcc_params <- paste0(dcc_D_params, '|', dcc_Q_params)

  ## Meta-data needed for printing
  ## TODO: Revisit this; some can be removed. Kitchen sink for now.
  metaNames <- c("param", "distribution", "num_dist", "iter", "chains",
                 "elapsed_time", "date", "nt", "TS_names", "mgarchQ",
                 "mgarchP", "meanstructure", "sampling_algorithm")
  meta <- with(object, mget(metaNames))
  meta$xC <- !all(object$xC == 0)
  out <- list()
  out$meta <- meta
  out$meta$digits <- digits
  
  ## Get the model summaries. print.summary.dcnet will process this + meta.
  params <- switch(object$param,
                   CCC = paste0(ccc_params,'|', var_params, '|', common_params),
                   DCC = paste0(dcc_params,'|',var_params,'|', common_params),
                   NULL
                   )
  if(is.null(params)) {
    stop("dcnet object 'param' does not match a supported model. ",
         object$param, "is not one in ", paste0(supported_models, collapse = ", "), ".")
  }


  out$model_summary <-  .get_stan_summary(object, params, CrI,
                                          sampling_algorithm = object$sampling_algorithm)
  
  class(out) <- "summary.bmgarch"
  return(out)
}

##' Internal function to compute summary stats from posterior draws
##' Called in .get_stan_summary
##' @title 
##' @param X 
##' @param CrI 
##' @return 
##' @author philippe
.summarize_draws <- function(X, CrI ) {
  mn <- mean(X)
  sd <- sd(X)
  md <- median(X )
  qf <- quantile(X, probs = CrI )
  return(c(mean =  mn,  sd =  sd, mdn =  md, qf) )
}



.get_stan_summary <- function(object, params, CrI, sampling_algorithm ) {
  if(sampling_algorithm == 'HMC') {
    cols <- c("mean","sd",paste0(CrI*100, "%"), "n_eff", "Rhat")            
  } else { ## variational does not return n_eff and Rhat
    if(sampling_algorithm == 'variational' ) {
      cols <- c("mean","sd",paste0(CrI*100, "%"))                
    }
  }

  ## Depending on estimation method, extraction differs:
  ## Variational
  if( object$sampling_algorithm == "variational" ) {

    ## extracting draws from here can take a long time -- awkward user experience?
    ## params contains variables of interest given the model type (eg., CCC or DCC)
    ## col_selct returns column position of given variable in params
    col_select <-
      which( grepl(params, colnames(object$model_fit$draws( ))) )
    
    ## subset draws of fitted model
    draws <-
      object$model_fit$draws( )[, col_select]

    ## Summarize draws
    model_summary <-
      t( apply(draws,  2,  FUN = .summarize_draws, CrI ) )

    return(model_summary)
      
  } else if(object$sampling_algorithm == "HMC" ) {
    print("Todo")
  } 
}


##' @title Print method for dcnet.summary object
##' @param x 
##' @param ... 
##' @return x (invisible) 
##' @author philippe
print.summary.dcnet <- function(x,  ... ) {
  if(x$meta$param == "CCC") {
    .print.summary.ccc(x)
  } else if(x$meta$param == "DCC") {
    .print.summary.dcc(x)
  } 
  .newline(2)

  .print.summary.means(x)
  .newline(2)

  ## if(x$meta$xC) {
  ##   .print.summary.beta(x)
  ##   .newline(2)
  ## }

  ## if(x$meta$num_dist == 1) {
  ##   .print.summary.nu(x)
  ##   .newline(2)
  ## }

  .print.summary.lp(x)
  .newline()

  return(invisible(x))
}

.print.summary.ccc <- function(bmsum, ...) {
   # Meta-data
    cat(paste0("Model: mlVAR-", bmsum$meta$param, "\n"))
    cat("Basic Specification: ")
    cat("H_t = D_t R D_t")
    .newline()
    .tab()
    cat("diag(D_t) = sqrt(h_[ii,t]) = c_h + a_h*y^2_[t-1] + b_h*h_[ii, t-1")
    .newline()

    # Sampling configuration
  
    .print.config(bmsum)
    # Get indices for params
    ms <- bmsum$model_summary
    ms <- ms[!grepl("c_h\\[", rownames(ms)),] # Remove c_h; will print c_h_var
    garch_h_index <- grep("_h", rownames(ms))
    cond_corr_index <- grep("R", rownames(ms))

    # Settings
    nt <- bmsum$meta$nt
    P <- bmsum$meta$mgarchP
    Q <- bmsum$meta$mgarchQ
    digits <- bmsum$meta$digits

    # Shortened TS names, if needed.
    short_names <- abbreviate(bmsum$meta$TS_names, minlength = 2)
    # Off-diagonals
    cormat_index <- matrix(1:(nt*nt), nrow = nt)
    corr_only <- cormat_index[lower.tri(cormat_index)]
    diag_only <- diag(cormat_index)
    ## obtain all combinations of TS varnames for A and B in BEKK
    full_varnames <- expand.grid( short_names, short_names)
    ## obtain off-diagonal TS varnames
    od_varnames <- full_varnames[corr_only, ]

    cat(paste0("GARCH(", P, ",", Q, ")"), " estimates for conditional variance:")
    .newline(2)


    #########
    # GARCH #
    #########
}



##' @title Print helper for DCC.
##' @param bmsum summary.bmgarch object.
##' @return Void.
##' @author Stephen R. Martin, Philippe Rast
##' @keywords internal
.print.summary.dcc <- function(bmsum) {
                                        # Meta-data
  cat(paste0("Model: mlVAR-", bmsum$meta$param, "\n"))
  cat("Basic Specification: ")
  cat("H_t = D_t R D_t")
  .newline()
  .tab()
  cat("diag(D_t) = sqrt(h_ii,t) = c_h + a_h*y^2_[t-1] + b_h*h_[ii,t-1]")
  .newline()
  .tab()
  cat("R_t = Q^[-1]_t Q_t Q^[-1]_t = ( 1 - a_q - b_q)S + a_q(u_[t-1]u'_[t-1]) + b_q(Q_[t-1])")
  .newline()

  ## Sampling configuration
  .print.config(bmsum)

  ## Get indices for params
  ms <- bmsum$model_summary
  ##    ms <-  ms[!grepl("c_h\\[", rownames(ms)),] # Remove c_h; will print c_h_var
  garch_h_index <- grep("_h", rownames(ms))
  garch_q_index  <- grep("_q", rownames(ms) )
  cond_corr_index <- grep("R", rownames(ms))
  S_index = grep("S", rownames(ms))

  ## Settings
  nt <- bmsum$meta$nt
  P <- bmsum$meta$mgarchP
  Q <- bmsum$meta$mgarchQ
  digits <- bmsum$meta$digits

  ## Shortened TS names, if needed.
  short_names <- abbreviate(bmsum$meta$TS_names, minlength = 2)
  ## Off-diagonals
  cormat_index <- matrix(1:(nt*nt), nrow = nt)
  corr_only <- cormat_index[lower.tri(cormat_index)]
  diag_only <- diag(cormat_index)
  ## obtain all combinations of TS varnames for A and B in BEKK
  full_varnames <- expand.grid( short_names, short_names)
  ## obtain off-diagonal TS varnames
  od_varnames <- full_varnames[corr_only, ]

#########
                                        # GARCH #
#########
  cat(paste0("GARCH(", P, ",", Q, ")"), " estimates for conditional variance on D:")
  .newline(2)

  rn <- rownames(ms[garch_h_index,])
  for ( i in 1:nt ) {
    replace <- grep(paste0( as.character(i), "\\]"), rn) 
    rn[replace] <- gsub(paste0( as.character(i), "\\]" ), paste0(short_names[i]), rn)[replace]
  }
  rn <- gsub("\\[", "_", rn)

  ## Save into new object to change rownames
  garch_h_out <- ms[garch_h_index,]
  ## Assign new rownames
  rownames(garch_h_out) <- rn
  ## print garch parameters
  print(round( garch_h_out, digits = digits) )
  .newline(2)

#####
                                        # Q #
#####
  cat("GARCH(1,1) estimates for conditional variance on Q:")
  .newline(2)
  rn = rownames(ms[garch_q_index,])
  for ( i in 1:nt ) {
    replace = grep(paste("\\[", "\\]", sep=as.character(i)), rn)
    rn[replace] = gsub(paste("\\[", "\\]", sep=as.character(i)), paste0("_", short_names[i]), rn)[replace]
  }

  ## Save into new object to change rownames
  garch_q_out = ms[garch_q_index,]
  ## Assign new rownames
  rownames(garch_q_out) = rn
  ## print garch parameters
  print(round( garch_q_out, digits = digits) )
  .newline(2)

  cat("Unconditional correlation 'S' in Q:")
  .newline(2)
  S_out <- ms[S_index[corr_only],]
  if (nt == 2) {
    tmp <- matrix( S_out, nrow = 1 )
    rownames(tmp) <- paste( paste0("R_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
    colnames(tmp) <- names(S_out)
    S_out <- tmp 
  } else {
    rownames(S_out) <- paste( paste0("S_", substr(od_varnames[ ,1], 1, 2) ), substr(od_varnames[ ,2], 1, 2) , sep = '-')
  }

  print(round(S_out, digits = digits))
}





##' @title Print helper for means component.
##' @param bmsum summary.bmgarch object.
##' @return Void.
##' @author Stephen R. Martin, Philippe Rast
##' @keywords internal
.print.summary.means <- function(bmsum) {
  ms <- bmsum$model_summary
  nt <- bmsum$meta$nt
  TS_names <- bmsum$meta$TS_names

    intercept_index <- grep("phi0", rownames(ms))
    ar_index <- grep("phi\\[", rownames(ms))
    ma_index <- grep("theta\\[", rownames(ms))

    ar_indices <- gsub("phi\\[([[:digit:]]+),([[:digit:]]+)]", "\\1,\\2", rownames(ms)[ar_index])
    ma_indices <- gsub("theta\\[([[:digit:]]+),([[:digit:]]+)]", "\\1,\\2", rownames(ms)[ma_index])
    ar_indices <- strsplit(ar_indices, ",")
    ma_indices <- strsplit(ma_indices, ",")
##    ar_indices <- apply(do.call(rbind, ar_indices), 1:2, as.numeric)
##    ma_indices <- apply(do.call(rbind, ma_indices), 1:2, as.numeric)

    rownames(ms)[intercept_index] <- paste0("(Intercept)_", bmsum$meta$TS_names)
##    rownames(ms)[ar_index] <- paste0("Phi_", TS_names[ar_indices[,1]], "-", TS_names[ar_indices[,2]])
##    rownames(ms)[ma_index] <- paste0("Theta_", TS_names[ma_indices[,1]], "-", TS_names[ma_indices[,2]])

    if(bmsum$meta$meanstructure == 0) {
        means_index <- intercept_index
        msg <- "Intercept estimates on the location:"
    } else if(bmsum$meta$meanstructure == 1) {
        means_index <- c(intercept_index, ar_index, ma_index)
        msg <- "ARMA(1,1) estimates on the location:"
    } else if(bmsum$meta$meanstructure == 2) {
        means_index <- c(intercept_index, ar_index)
        msg <- "VAR(1) estimates on the location:"       
    }
    cat(msg)
    .newline(2)
    print(round(ms[means_index,], digits = bmsum$meta$digits))
}


##' @title Print helper for Sampling Config.
##' @param bmsum summary.bmgarch object.
##' @return Void.
##' @author Philippe Rast
##' @keywords internal
.print.config <- function(bmsum) {
  .newline()
  cat("Sampling Algorithm: ", bmsum$meta$sampling_algorithm)
  .newline()
  cat("Distribution: ", bmsum$meta$distribution)
  .newline()
  .sep()
  cat("Iterations: ", bmsum$meta$iter)
  .newline()
  cat("Chains: ", bmsum$meta$chains)
  .newline()
  cat("Date: ", bmsum$meta$date)
  .newline()
  cat("Elapsed time (min): ", round((max(bmsum$meta$elapsed_time))/60, 2))
  .newline(2)
}


##' @title Print helper - Separator, new line
##' @return Prints "---" and a new line.
##' @keywords internal
.sep <- function() {
  cat("---")
  .newline()
}

##' @title Print helper - tab
##' @return Prints tab.
##' @keywords internal
.tab <- function() {
  cat("\t")
}

##' @title Print helper - Return new line(s).
##' @param n Integer (Default: 1). Number of new lines.
##' @return Prints new lines.
##' @keywords internal
.newline <- function(n = 1) {
  for(i in 1:n) {
    cat("\n")
  }
}

##' @title Print helper for LP component.
##' @param bmsum summary.bmgarch object.
##' @return Void.
##' @keywords internal
.print.summary.lp <- function(bmsum) {
  cat("Log density posterior estimate:")
  .newline(2)
  print(round(bmsum$model_summary["lp__",], digits = bmsum$meta$digits))
}
