# Check log-concave

is_logconcave <- function(d,x_lo,x_hi,twice_differentiable=TRUE) {
  is_logc<-TRUE
  d_log <- function(x) log(d(x))
  
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
