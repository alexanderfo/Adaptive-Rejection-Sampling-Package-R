# Find the mode of the density within boundaries

find_mode <- function(density, lb = -Inf, ub = Inf) {
  if (lb != -Inf & ub != Inf) {
    out <- optim((lb+ub)/2, density, method="Brent", lower=lb, upper=ub, 
                 control = list(fnscale=-1))
  } else {
    out <- optim(0, density, method="BFGS", control = list(fnscale=-1))
    if (out$convergence == 1) stop("Density not log-concave: no maxima")
  }
  mode <- out$par
  h_mode <- out$value
  return(c(mode, h_mode))
}