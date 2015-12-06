# Initialization, when only upper and lower bound for x is provided and log(density)=log(f(x))=h(x)

setwd("~/src/stat243-project/ars/R")
source("evaluate_deriv.R")
install.packages("micEcon")
library(micEcon) #
?translogCheckCurvature


intialize<-function(h,x_lo,x_hi){
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
  <<
  colnames(vertices)<-c("x","h(x)","h_prime(x)","secant")
  rownames(vertices)<-NULL # remove row names created by rbind
  
  vertices[1,4]<-calc_secant(vertices,1,2)
  
  
  # Initilize u
  
  # Can also be used to update u. Returns the function for u (the tangent), for x at index
  create_u <- function(vertices,index,h) {
    u<- function (x) vertices[index,2]+(x-vertices[index,1])*evaluate_deriv(h,vertices[index,1])
    return(u)
  }
  
  u1<-create_u(vertices,1)
  u2<-create_u(vertices,2)
  
  # Exponentiate a function, e.g. the u function created aboves, returns an exponentiated function 
  exp_fun <- function(fun){
    new_exp_fun<- function (x) exp(fun(x))
    return(new_exp_fun)
  }
  
  
  # Initialize z, only creates first value of z, can be made into a function for all other z's
  z1<-(h(x_hi)-h(x_lo)-x_hi*evaluate_deriv(h,x_hi)+x_lo*evaluate_deriv(h,x_lo))/(evaluate_deriv(h,x_lo)-evaluate_deriv(h,x_hi))
  
  
  # Initialize l
  
  #Can also be used to update l. Returns the function for l (the tangent), for x at index
  # NB. for x<x_lo and x>x_hi then l(x) = +-Inf
  # The secant function between x_lo to x_hi
  create_l <- function(vertices,index) {
    l<-function(x) vertices[index,4]*(x-vertices[index,1])+vertices[index,2]
    return(l)
  }
  
  
  
  
  
}