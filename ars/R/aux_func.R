# Create the l functions (secant lines functions connecting the x's)
# Returns the function for l, for x at index
# NB. for x<x_lo and x>x_hi then l(x) = +-Inf
# The secant function between x_lo to x_hi
create_l <- function(vertices,index) {
  l<-function(x) vertices[index,4]*(x-vertices[index,1])+vertices[index,2]
  return(l)
}

# Create the u functions (upper hull)
# Returns the function for u (the tangent), for x at index
create_u <- function(vertices,index,h) {
  u<- function (x) vertices[index,2]+(x-vertices[index,1])*evaluate_deriv(h,vertices[index,1])
  return(u)
}

# Exponentiate a function, e.g. the u function created aboves, returns an exponentiated function 
exp_fun <- function(fun){
  new_exp_fun<- function (x) exp(fun(x))
  return(new_exp_fun)
}

# calculate the slope value of a secant line given the vertices matrix and the row indices 
calc_secant <- function(vertices, row_id1, row_id2){
  numerator <- vertices[row_id2,2] - vertices[row_id1, 2]
  denominator <- vertices[row_id2,1] - vertices[row_id1, 1]
  return(numerator/denominator)
}

# given an arbitrary x value, find the index(rank) of the next bigger z value in func_list
# may be extended to find not just z, but also x
find_bin <- function(x, hi_vec){
  row <- sapply(1:length(x), function(i) which(x[i] <= hi_vec)[1])
  return(row)
}