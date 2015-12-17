#' @title Check convergence of a density
#' 
#' @description Check if a given density has finite integral over the sampling range. If it diverges, error is reported.
#'
#' @param density The unnormalized target density
#' @param lower/upper The lower / upper bound of the domain, where the points are sampled from
#' 
#' @return A finite number the unnormalized density converge to.
check_density_convergence <- function(density, lower, upper) {  
  integral <- try(integrate(density, lower, upper), silent = TRUE)
  if (is(integral, "try-error")) {
    stop("Bad density: density diverges")
  }
  return(integral$value)
}