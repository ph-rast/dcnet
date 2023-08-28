get_target_stan_path <- function() {
  package_path <-
    system.file(package = "dcnet", "stan")
  return(package_path)
}
