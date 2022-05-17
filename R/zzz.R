.onLoad <- function(libname, pkgname) {
  ## Check for presence of cmdstanr
  if( !require("cmdstanr") ) stop("The `cmdstanr` library as well as CmdStan are required to run this pacakge. Check out https://mc-stan.org/cmdstanr")
}
