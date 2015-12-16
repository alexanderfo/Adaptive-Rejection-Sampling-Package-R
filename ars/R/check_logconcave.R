source("ars/R/evaluate_deriv.R")

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

check_local_concave <- function(u, l, xvec, halfrange = 1e-16){
  if(any((u(xvec+halfrange)) < l(xvec+halfrange))) return(FALSE)
  if(any((u(xvec-halfrange)) < l(xvec-halfrange))) return(FALSE)
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