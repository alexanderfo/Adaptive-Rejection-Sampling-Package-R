# Create the l functions (secant lines)
# Returns the function for l (the tangent), for x at index
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

# calculate the secant given the vertices matrix and the row indices 
calc_secant <- function(vertices, row_id1, row_id2){
  numerator <- vertices[row_id2,2] - vertices[row_id1, 2]
  denominator <- vertices[row_id2,1] - vertices[row_id1, 1]
  return(numerator/denominator)
}

# log-concavity ensures that the denominator is never 0 analytically
calc_intersection <- function(h, vertices, idx1, idx2){
  x_lo <- vertices[idx1,1]
  x_hi <- vertices[idx2,1]
  numerator <- h(x_hi) - h(x_lo) - x_hi * evaluate_deriv(h, x_hi) + x_lo * evaluate_deriv(h, x_lo)
  denominator <- evaluate_deriv(h, x_lo) - evaluate_deriv(h, x_hi)
  return(numerator/denominator)
}

# Mengfei & Daisy
draw_sample <- function(envelope, z, num_of_samples = 1){
	# Draw sampling point candidates
  # Args: 
  #   envelope: vector of exponentiated linear envelope functions
  #   z: vector of intersection points
  # Output:
  #   candidates: vector sampling candidates
  numerator <- sapply(1:length(envelope), function(i) integrate(envelope[[i]], z[i], z[i+1])$value)
  denominator <- sum(numerator)
  cdf_end_pts <- cumsum(numerator/denominator)
  w <- runif(num_of_samples)
  # locate the bin each w falls in, i.e. find i
  w_idx <- sapply(w, function(i) which(i <= cdf_end_pts)[1])
  f <- function(x, i, w_i) {
    # Args:
    #   x: one sampling candidate
    #   i: the index of the Uniform variable in the w vector
    #   w_i: the bin index that the Uniform variable falls in
    integrate(envelope[[w_i]], lower = z[w_i], upper = x)$value/denominator - w[i] + ifelse(w_i == 1, 0, cdf_end_pts[w_i-1])
    }
  candidates <- sapply(1:num_of_samples, 
                       function(i) {
                         uniroot(f,interval = c(z[w_idx[i]], z[w_idx[i]+1]), 
                                 i, w_idx[i], tol = 1e-6)$root})
  return(candidates)
}

find_bin <- function(x, z_hi_vec){
  cand <- sapply(1:length(x), function(i) which(x[i] <= z_hi_vec)[1])
  return(cand)
}

# little test
f1 <- function(x) dnorm(x)
f2 <- function(x) dexp(x)
envelope <- c(f1, f2)
z <- c(-1, 0.5, 2)
# see code efficiency
tmp <- tempfile()
Rprof(tmp, interval = 0.1)
a=draw_sample(envelope, z, num_of_samples = 5000)
Rprof(NULL)
summaryRprof(tmp)

