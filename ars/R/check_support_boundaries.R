#' Check the boundary is supportive to the function f
#' @param f The density function (which could be unnormalized)
#' @param lower/upper The boundaries provided by the user
#' @return If the unnormalized density is positive and finite within the bounds; otherwise, the program is aborted
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
  return(TRUE)
}