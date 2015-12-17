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

# visually display the log-density, the upper hull, and lower squeezing function
plot_iter <- function(h, u, l, lb, ub){
  if(is.infinite(lb)) lb <- -38
  if(is.infinite(ub)) ub <- 38
  dat <- seq(lb, ub, by=0.1)
  plot(dat, h(dat), type='l')
  lines(dat, u(dat), lty=5)
  lines(dat, l(dat), lty=3)
}

# check secant values are all correctly calculated in the vertices matrix
check_secant <- function(vertices){
  eps <- 1e-8
  len <- length(vertices[,1])
  for(i in 1:(len-1)){
    if(vertices[i,4] * (vertices[i+1,1] - vertices[i,1]) - (vertices[i+1,2] - vertices[i, 2]) > eps){
      stop(paste("The secant value is incorrect in ", i, "-th row."))
    }
  }
  print("All pass!")
}

# visually compare the sampled density with the true density
compare_densities<-function(x_real,x_ars) {
  print(density(x_real))
  print(density(x_ars))
  
  par(mfrow=c(1,1))
  plot(density(x_real),type="l",col="green",main="real dens green, approx dens red",xlab="x")
  lines(density(x_ars),type="l",col="red")
}

v_integrate <- Vectorize(integrate, vectorize.args = c("lower", "upper"))