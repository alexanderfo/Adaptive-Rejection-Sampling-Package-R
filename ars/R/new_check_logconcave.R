source("ars/R/evaluate_deriv.R")

is_logconcave <- function(h, x_lo, x_hi, mode, ...) {
  # h is the density function defined elsewhere
  # ...: arguments to be passed to h
  h_x_lo <- h(x_lo, ...)
  h_x_hi <- h(x_hi, ...)
  h_mode <- h(mode, ...)
  
  eps <- sqrt(.Machine$double.eps)
  
  # a pre-check on the input density and boundaries
  # 1 - uniform; 2 - mode is equal to lb; 3 - mode is equal to ub; 4 - lb, ub, and mode are spread afar
  if ((abs(h_mode - h_x_lo) < eps) && (abs(h_x_hi - h_mode) < eps)) return(1)
  else if (abs(mode - x_lo) < eps) return(2)
  else if (abs(mode - x_hi) < eps) return(3)
  else return(4)
}

is_logconcave_core <- function(h,x_lo,x_hi,twice_differentiable=TRUE, ...) {
  is_logc<-TRUE
  d_log <- function(x, ...) h(x, ...)
  
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
# 
# is_logconcave_core <- function(d, start, end) {
#   x_vals <- seq(start, end, length.out = 1000)
#   vals <- sapply(x_vals, d)
#   slopes <- try(diff(vals) / (x_vals[2] - x_vals[1]))
#   print(slopes)
#   if (all(diff(slopes) <= 0)) {
#     return(TRUE)
#   }
#   return(FALSE)
# }
# 
# d <- dnorm
# h <- function(x) log(d(x))
# is_logconcave(h, x_lo = -5, x_hi = 5, mode = 0, h_mode = h(0))
# is_logconcave_core(h, -5, 5, TRUE)
# d <- dunif
# h <- function(x, ...) log(d(x, ...))
# is_logconcave(h, x_lo = -5, x_hi = 5, mode = 0, h_mode = 1, -5, 5)
# undebug(is_logconcave)
# d <- dexp
# h <- function(x, ...) log(d(x, ...))
# is_logconcave(h, 0, 5, mode = 0, 2)
