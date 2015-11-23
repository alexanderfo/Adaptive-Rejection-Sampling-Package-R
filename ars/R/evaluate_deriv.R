## Function that numerically evaluates derivative of a function / density


evaluate_deriv <- function(h,x,diff=10^-6) {
  x_diff=c(x-diff,x+diff)
  
  der_value<-diff(h(x_diff))/diff(x_diff)
  
  return(der_value)
}
