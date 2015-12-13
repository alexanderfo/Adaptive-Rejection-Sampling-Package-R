# Find the mode of the density within boundaries
# density is in log scale
find_mode <- function(density, lb = -Inf, ub = Inf) {
  if (!is.infinite(density(lb)) && !is.infinite(density(ub)) && !is.infinite(lb) && !is.infinite(ub)) {
    out <- optim((lb+ub)/2, density, method="Brent", lower=lb, upper=ub, 
                 control = list(fnscale=-1))
  } else {
    out <- optim(0, density, method="BFGS", control = list(fnscale=-1))
    if (out$convergence == 1) stop("Density not log-concave: no maxima")
  }
  mode <- out$par
  f_mode <- out$value
  return(c(mode, f_mode))
}