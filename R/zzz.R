.onLoad <- function(libname, pkgname) {
  ## Check for presence of cmdstanr
  if( !require("cmdstanr") ) stop("The `cmdstanr` library as well as CmdStan are required to run this pacakge. Check out https://mc-stan.org/cmdstanr")

  message("Compiling Stan models...")
  
  ## Set the path to your models
  stan_path <- system.file("stan", package = pkgname)

  ## Paths to the Stan model files
  ccc_file <- file.path(stan_path, "VAR.stan")
  dcc_file <- file.path(stan_path, "DCCMGARCHrandQ.stan")
  dccr_file <- file.path(stan_path, "DCCMGARCHrandS.stan")

  ## Compile the models if not already compiled
  assign("ccc_model", cmdstanr::cmdstan_model(ccc_file), envir = .GlobalEnv)
  assign("dcc_model", cmdstanr::cmdstan_model(dcc_file), envir = .GlobalEnv)
  assign("dccr_model", cmdstanr::cmdstan_model(dccr_file), envir = .GlobalEnv)
  
  message("Stan models compiled successfully.")

}
