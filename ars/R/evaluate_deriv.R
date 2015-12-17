#' @title Calculates first order derivative
#' 
#' @description Numerically calculate the derivative locally for a continuous differentiable function
#'
#' @param h A differentiable function
#' @param x The abscissa at which the derivative is taken
#' @param diff The infinidesimal piece for calculation
#' 
#' @return The derivative value at point x for function h
evaluate_deriv <- function(h, x, diff=10^-12) {
  x_diff <- c(x-diff, x+diff)
  der_value <- diff(h(x_diff)) / diff(x_diff)
  return(der_value)
}

# Source: http://blog.quantitations.com/tutorial/2013/02/12/numerical-derivatives-in-r/

#' @title Calculates the derivatives at specified order
#' 
#' @description Calculate the derivative locally for a function of specific order
#'
#' @param f The function to take derivative for
#' @param x The abscissa at which the derivative is taken
#' @param order The order of the derivative
#' @param delta The infinidesimal piece for calculation
#' 
#' @return The derivative value at point x for function h of the specified order
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