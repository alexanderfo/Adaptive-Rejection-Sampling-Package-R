# check that the secants are calculated correctly
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

# Test for der-function
source("evaluate_deriv.R")
f<-function(x) return(x^3) # fun
evaluate_deriv(f,x=3) #OK