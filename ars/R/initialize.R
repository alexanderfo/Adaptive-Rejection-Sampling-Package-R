# Initialization, when only upper and lower bound for x is provided and log(density)=log(f(x))=h(x)

setwd("~/git/stat243-project/ars/R")
source("evaluate_deriv.R")
install.packages("micEcon")
library(micEcon) #
?translogCheckCurvature

# Initialize l
# Can also be used to update l. Returns the function for l (the tangent), for x at index
# NB. for x<x_lo and x>x_hi then l(x) = +-Inf
# The secant function between x_lo to x_hi
create_l <- function(vertices,index) {
  l<-function(x) vertices[index,4]*(x-vertices[index,1])+vertices[index,2]
  return(l)
}

# Initilize u
# Can also be used to update u. Returns the function for u (the tangent), for x at index
create_u <- function(vertices,index,h) {
  u<- function (x) vertices[index,2]+(x-vertices[index,1])*evaluate_deriv(h,vertices[index,1])
  return(u)
}

# Exponentiate a function, e.g. the u function created aboves, returns an exponentiated function 
exp_fun <- function(fun){
  new_exp_fun<- function (x) exp(fun(x))
  return(new_exp_fun)
}

calc_secant <- function(vertices, row_id1, row_id2){
  numerator <- vertices[row_id2,2] - vertices[row_id1, 2]
  denominator <- vertices[row_id2,1] - vertices[row_id1, 1]
  return(numerator/denominator)
}

calc_intersection <- function(h, vertices, idx1, idx2){
  x_lo <- vertices[idx1,1]
  x_hi <- vertices[idx2,1]
  numerator <- h(x_hi) - h(x_lo) - x_hi * evaluate_deriv(h, x_hi) + x_lo * evaluate_deriv(h, x_lo)
  denominator <- evaluate_deriv(h, x_lo) - evaluate_deriv(h, x_hi)
  return(numerator/denominator)
}

init_vertices<-function(h,x_lo,x_hi){
  # defensive programming, part check log-concavity
  if(is.infinite(x_lo)) {
    x_lo <- -10^16
    i_lo=16 #exponent for -10^i_lo
    while(evaluate_deriv(h,x_lo)<=0) {
      x_lo <- x_lo - 10^i_lo
      i_lo<-i_lo+1
      if(is.infinite(x_lo)) stop("NOT LOG-CONCAVE")
    }
  }
  
  if(is.infinite(x_hi)) {
    x_hi <- 10^16
    i_hi=16
    while(evaluate_deriv(h,x_lo)>=0) {
      x_hi <- x_hi + 10^i_hi
      i_hi<-i_hi+1
      if(is.infinite(x_hi)) stop("NOT LOG-CONCAVE")
    }
  }
  
  # Build matrix vertices that defines: x values, h(x) values, h_prime(x) value, and the secant slope between x1 to x2 (stored at index 1))
  row1<-c(x_lo,h(x_lo),evaluate_deriv(h,x_lo),NA)
  row2<-c(x_hi,h(x_hi),evaluate_deriv(h,x_hi),NA)
  
  vertices<-rbind(row1,row2)
  colnames(vertices)<-c("x","h(x)","h_prime(x)","secant")
  rownames(vertices)<-NULL # remove row names created by rbind
  
  vertices[1,4]<-calc_secant(vertices,1,2)
  
  return(vertices)
}

# input_lo & input_hi are the user-input range of the density
# these are used as z_0 and z_k for initalization
init_piecewise <- function(vertices, h, input_lo, input_hi){
  u1 <- create_u(vertices, 1, h)
  u2 <- create_u(vertices, 2, h)
  l1 <- create_l(vertices, 1)
  
  u_list <- c(u1,u2)
  exp_u_list <- list(exp_fun(u1), exp_fun(u2))
  l_list <- c(l1)
  exp_l_list <- list(exp_fun(l1))
  
  # Initialize z, only creates first value of z, can be made into a function for all other z's
  z1 < -calc_intersection(h, vertices, 1, 2)
  z_lo_vec <- c(input_lo, z1)
  z_hi_vec <- c(z1, input_hi)
  
  x_lo_vec <- c(vertices[1,1],NA)
  x_hi_vec <- c(vertices[2,1],NA)
  
  func_list <- list(u=c(u1,u2), exp_u=exp_u_list, z_lo=z_lo_vec, z_hi=z_hi_vec,
                     l=l_list, exp_l=exp_l_list, x_lo=x_lo_vec, x_hi=x_hi_vec)
  return(func_list)
}
