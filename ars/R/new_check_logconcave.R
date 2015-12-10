source("ars/R/evaluate_deriv.R")

is_logconcave <- function(h, x_lo, x_hi, mode, h_mode) {
  # h is the density function defined elsewhere
  h_x_lo <- h(x_lo)
  h_x_hi <- h(x_hi)
  
  if (((h_mode - h_x_lo) < .Machine$double.eps) && ((f_x_hi - f_mode) < .Machine$double.eps)) {
    # uniform distribution
    warning("Uniform distribution: runif is used to generate sample")
    return(1)
  } else if (mode - x_lo < .Machine$double.eps) {
    # exponential like distribution: mode = x_lo
    if (is_logconcave_core(h, mode, x_hi)) {return(2)}
  } else {
    # x_lo < mode < x_hi
    left_logc <- is_logconcave_core(h, x_lo, mode, TRUE)
    right_logc <- is_logconcave_core(h, mode, x_hi, TRUE)
    if (left_logc && right_logc) {return(3)}
  }
  return(FALSE)
}

is_logconcave_core <- function(h,x_lo,x_hi,twice_differentiable=TRUE) {
  is_logc<-TRUE
  d_log <- h
  
  if(twice_differentiable) {
    steps<-(x_hi-x_lo)/1000
    x_vals<-round(seq(x_lo,x_hi,steps),5)
    
    #f_vals<-sapply(x_vals, d_log) #check function values
    #if( !all( !is.nan(f_vals) ) ) {
    #is_logc<-FALSE
    #return(is_logc)
    #}
    
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


############################################
'''
is_logconcave_core <- function(d, start, end) {
  x_vals <- seq(start, end, length.out = 1000)
  vals <- sapply(x_vals, d)
  slopes <- try(diff(vals) / (x_vals[2] - x_vals[1]))
  print(slopes)
  if (all(diff(slopes) <= 0)) {
    return(TRUE)
  }
  return(FALSE)
}

d <- dnorm
h <- function(x) log(d(x))
is_logconcave(h, x_lo = -5, x_hi = 5, mode = 0, h_mode = h(0))
is_logconcave_core(h, -5, 5, TRUE)
'''