source("ars/R/find_mode.R")

check_support_boundaries <- function(f, lower, upper) {
  # Args:
  #   f: unnormalized target density
  #   lower: lower bound of the unnormalized target density
  #   upper: upper bound of the unnormalized target density
  # Returns:
  #   TRUE if the unnormalized density is positive within bounds
  
  # check lower smaller than upper
  if (lower >= upper) {stop("Lower bound greater than or equal to upper bound")}
  
  # check positiveness of unnormalized density at bounds
  f_lval <- f(lower)
  f_uval <- f(upper)
  if (is.infinite(f_lval) || f_lval < 0 || is.infinite(f_uval) || f_uval < 0) {
    stop("Bad bounds: density is +/-Inf or negative at upper/lower bounds")
  }
  
  # check the bounds enclose the mode
  mode <- find_mode(f, lb = lower, ub = upper)
  if (f_lval < .Machine$double.eps && f_uval < .Machine$double.eps &&
      mode[2] - f_lval < .Machine$double.eps) {
    stop("Bad bounds: density is 0 everywhere within bounds")
  }
  
  # check boundness of the bounds on the density
  #if (f_lval == 0) {
  #  stop("Lower bound does not bound the density correctly")
  #} else if (f_uval == 0) {
  #  stop("Upper bound does not bound the density correctly")
  #}
  return(TRUE)
}