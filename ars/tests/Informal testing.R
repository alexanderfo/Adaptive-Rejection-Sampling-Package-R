## Informal testing

setwd("~/src/stat243-project")

# Test for der-function
source("evaluate_deriv.R")


f<-function(x) return(x^3) # fun

evaluate_deriv(f,x=3) #OK


# Test update_matrix-function

a<-matrix(rnorm(20),5,4)

update_matrix(a,0,rnorm(1),rnorm(1)) #OK
