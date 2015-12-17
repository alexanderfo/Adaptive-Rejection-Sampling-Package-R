#' @title Check the sampling condition based on shape
#' 
#' @description Decide the sampling condition of the procedure based on the shape and implied log-concavity
#'
#' @param h The log-density
#' @param x_lo/x_hi The lower and upper bounds of the sampling range
#' @param mode The abscissa of the mode of the density in the given range
#' 
#' @return The sampling condition coded in number 1-4. 1 - uniform distribution, 2 - left truncated, 3 - right truncated, and 4 - a regular density
is_logconcave_shape <- function(h, x_lo, x_hi, mode, ...) {
  # h is the density function defined elsewhere
  # ...: arguments to be passed to h
  h_x_lo <- h(x_lo, ...)
  h_x_hi <- h(x_hi, ...)
  h_mode <- h(mode, ...)

  eps <- 1e-6

  # a pre-check on the input density and boundaries
  # 1 - uniform; 2 - mode is equal to lb; 3 - mode is equal to ub; 4 - lb, ub, and mode are spread afar
  if ((abs(h_mode - h_x_lo) < eps) && (abs(h_x_hi - h_mode) < eps)) return(1)
  else if (abs(mode - x_lo) < eps) return(2)
  else if (abs(mode - x_hi) < eps) return(3)
  else return(4)
}

#' @title Check local concavity
#' 
#' @description Check the local concavity of a function by taking left/right points of a vector of chosen abscissae, and compare the tangent/secant values. If any of the inspected points has larger secant values than tagent, then the function is not concave.
#'
#' @param u The piecewise linear upper hull functions (tangent lines)
#' @param l The piecewise linear lower squeezing functions (secant lines)
#' @param xvec The abscissae of the chosen x values, the left/right of which are to be insepected
#' @param halfrange The distance of the inspected points from the chosen abscissae
#' 
#' @return FALSE if any tangent value is lower; TRUE if all tangent values are higher.

check_local_concave <- function(u, l, xvec, halfrange = 1e-6){
  eps <- 1e-8
  if(any(l(xvec+halfrange) - u(xvec+halfrange) > eps)) return(FALSE)
  if(any(l(xvec-halfrange) - u(xvec-halfrange) > eps)) return(FALSE)
  else return(TRUE)
}

is_logconcave_core <- function(h,x_lo,x_hi,twice_differentiable=TRUE, ...) {
  is_logc<-TRUE
  d_log <- function(x, ...) h(x, ...)

  if(twice_differentiable) {
    steps<-(x_hi-x_lo)/1000
    x_vals<-round(seq(x_lo,x_hi,steps),5)

    vals<-sapply(x_vals, der_fun, f = d_log, order = 2)

    if( !all( !is.nan(vals) ) ) {
      is_logc<-FALSE
      return(is_logc)
    }
    if(all(round(vals,6)<=0)) {
      return(is_logc)
    }
    else {
      is_logc<-FALSE
      return(is_logc)
    }
  }
}
