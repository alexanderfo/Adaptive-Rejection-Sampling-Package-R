#' Find the mode of the density within boundaries
find_mode <- function(density, lb = -Inf, ub = Inf) {
  if (!is.infinite(density(lb)) && !is.infinite(density(ub)) && !is.infinite(lb) && !is.infinite(ub)) {
    vals <- seq(lb, ub, by = .1)
    start_pt <- vals[which.max(density(vals))]
    out <- optim(start_pt, density, method="Brent", lower = max(start_pt-20, lb), upper = min(start_pt+20, ub), 
                 control = list(fnscale=-1))
  } else if (is.infinite(lb) && is.finite(ub)) {
    lb <- optim(min(-1,ub-1), function(x) {density(x) - 1e-18}, method = "BFGS")$par
    out <- optim(ub, density, method="Brent", lower=lb, upper=ub, 
                 control = list(fnscale=-1))
    lb <- optim(out$par-1, function(x) {density(x) - 1e-18}, method = "BFGS")$par
  } else if (is.finite(lb) && is.infinite(ub)) {
    ub <- optim(max(1,lb+1), function(x) {density(x) - 1e-18}, method = "BFGS")$par
    out <- optim(lb, density, method="Brent", lower=lb, upper=ub, 
                 control = list(fnscale=-1))
    ub <- optim(out$par+1, function(x) {density(x) - 1e-18}, method = "BFGS")$par
  } else {
    out <- optim(0, density, method="BFGS", control = list(fnscale=-1))
    if (out$convergence == 1) stop("Density not log-concave: no maxima")
    lb <- optim(out$par-1, function(x) {density(x) - 1e-18}, method = "BFGS")$par
    ub <- optim(out$par+1, function(x) {density(x) - 1e-18}, method = "BFGS")$par
  }
  mode <- out$par
  f_mode <- out$value
  return(c(mode, f_mode, lb, ub))
}