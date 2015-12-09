## Function that numerically evaluates derivative of a function / density
evaluate_deriv <- function(h, x, diff=10^-6) {
  x_diff <- c(x-diff, x+diff)
  der_value <- diff(h(x_diff)) / diff(x_diff)
  return(der_value)
}

# Source: http://blog.quantitations.com/tutorial/2013/02/12/numerical-derivatives-in-r/
der_fun <- function(f, x, ..., order = 1, delta = 10^-2, sig = 6) {
  # Numerically computes the specified order derivative of f at x
  vals <- matrix(NA, nrow = order + 1, ncol = order + 1)
  grid <- seq(x - delta/2, x + delta/2, length.out = order + 1)
  vals[1, ] <- sapply(grid, f, ...) - f(x, ...)
  for (i in 2:(order + 1)) {
    for (j in 1:(order - i + 2)) {
      stepsize <- grid[i + j - 1] - grid[i + j - 2]
      vals[i, j] <- (vals[i - 1, j + 1] - vals[i - 1, j])/stepsize
    }
  }
  return(signif(vals[order + 1, 1], sig))
}