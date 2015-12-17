# check integral of density < Inf
check_density_convergence <- function(density, lower, upper) {
  # Args:
  #   f: unnormalized target density
  #   lower: lower bound of the unnormalized target density
  #   upper: upper bound of the unnormalized target density
  # Returns:
  #    a finite number the unnormalized density converge to
  
  integral <- try(integrate(density, lower, upper), silent = TRUE)
  if (is(integral, "try-error")) {
    stop("Bad density: density diverges")
  }
  return(integral$value)
}