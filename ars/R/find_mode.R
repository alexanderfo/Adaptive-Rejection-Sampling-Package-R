#' @title Find mode of a density function
#' 
#' @description Calculate the abscissa of the mode of the input density function within the given range. It utilizes Brent when either one of the lower/upper bounds is finite, and BFGS when both of them infinite. The function also shrinks sampling range by assigning a finite value to the infinite bound(s) if there exist one.
#'
#' @param density The density function to sample from
#' @param lb/ub The range of the density function within which the mode is being seeked. The default values are set to plus/minus infinity.
#' 
#' @return A numeric vector of information found in the procedure. The 1st element is the abscissa of the mode in the range. The 2nd element is the density at the mode. The 3rd/4th elements are the updated(shrinked) lower/upper bound of the sampling range, if the input range is infinite.
find_mode <- function(density, lb = -Inf, ub = Inf) {

  # if both lb and ub are finite
  if (is.finite(lb) && is.finite(ub)) {
    vals <- seq(lb, ub, by = .1)
    start_pt <- vals[which.max(density(vals))]
    out <- optim(start_pt, density, method="Brent", lower = max(start_pt-20, lb), upper = min(start_pt+20, ub), 
                 control = list(fnscale=-1))
  } 

  # if lb is infinite and ub is finite, update lb accordingly
  else if (is.infinite(lb) && is.finite(ub)) {
    lb <- optim(min(-1,ub-1), function(x) {density(x) - 1e-18}, method = "BFGS")$par
    out <- optim(ub, density, method="Brent", lower=lb, upper=ub, 
                 control = list(fnscale=-1))
    lb <- optim(out$par-1, function(x) {density(x) - 1e-18}, method = "BFGS")$par
  } 

  # if lb is finite and ub is infinite, update ub accordingly
  else if (is.finite(lb) && is.infinite(ub)) {
    ub <- optim(max(1,lb+1), function(x) {density(x) - 1e-18}, method = "BFGS")$par
    out <- optim(lb, density, method="Brent", lower=lb, upper=ub, 
                 control = list(fnscale=-1))
    ub <- optim(out$par+1, function(x) {density(x) - 1e-18}, method = "BFGS")$par
  } 

  # if both lb and ub are infinite, update them both
  else {
    out <- optim(0, density, method="BFGS", control = list(fnscale=-1))
    if (out$convergence == 1) stop("Density not log-concave: no maxima")
    lb <- optim(out$par-1, function(x) {density(x) - 1e-18}, method = "BFGS")$par
    ub <- optim(out$par+1, function(x) {density(x) - 1e-18}, method = "BFGS")$par
  }

  # extract the mode and density at mode
  mode <- out$par
  f_mode <- out$value
  return(c(mode, f_mode, lb, ub))
}