# Test the ars
# setwd("~/git/stat243-project/")
source("ars/R/ars.R")
library(testthat)

compare_densities<-function(x_real,x_ars) {
  print(density(x_real))
  print(density(x_ars))
  
  par(mfrow=c(1,1))
  plot(density(x_real),type="l",col="green",main="real dens green, approx dens red",xlab="x")
  lines(density(x_ars),type="l",col="red")
}

n<-1000 #number of samples
################################################################################
############################ Log Concave Dist ##################################
#### Normal
# Case 1: Correct lower and upper bounds
test_that("ars correctly samples from standard normal", {
  x_real<-rnorm(n)
  x_ars<-ars(dnorm,n,-100,100)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  # compare_densities(x_real,x_ars) 
  expect_that(test$p.value >= 0.05, is_true())
})

# Case 2: No bounds
test_that("ars correctly samples from standard normal", {
  x_real<-rnorm(n)
  x_ars<-ars(dnorm,n)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  #compare_densities(x_real,x_ars) 
  expect_that(test$p.value >= 0.05, is_true())
})


# Case 3: lower bound is larger than upper bound
test_that("ars correctly samples from standard normal", {
  x_real<-rnorm(n)
  x_ars<-ars(dnorm,n,10,-10)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  #compare_densities(x_real,x_ars) 
  expect_that(test$p.value >= 0.05, is_true())
})
# Case 4: Can you check this?
#"Truncated distribution: the leftmost point is the mode."
#Error: Test failed: 'ars correctly samples from standard normal'
test_that("ars correctly samples from standard normal", {
  x_real<-rnorm(n)
  x_ars<-ars(dnorm,n,1,2)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  #compare_densities(x_real,x_ars) 
  expect_that(test$p.value >= 0.05, is_true())
})


################################################################################
##### Truncated normal

############################### Informal test ##################################

# library(truncnorm)
# x_real<-rtruncnorm(n,a=-10,b=10)
# x_ars<-ars(dnorm,n,-10,10)
# 
# compare_densities(x_real,x_ars) #OK
# test <- ks.test(x_ars, x_real) # p-value = 0.8127
# expect_that(test$p.value >= 0.05, is_true())

###############################################################################
# Case 1: Correct lower and upper bounds
test_that("ars correctly samples from truncated normal", {
  library(truncnorm)
  x_real<-rtruncnorm(n,-10,10)
  x_ars<-ars(dnorm,n,-10,10)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  #compare_densities(x_real,x_ars) 
  expect_that(test$p.value >= 0.05, is_true())
})

# Case 2: No bounds
test_that("ars correctly samples from truncated normal", {
  library(truncnorm)
  x_real<-rtruncnorm(n,-10,10)
  x_ars<-ars(dnorm,n)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  #compare_densities(x_real,x_ars) 
  expect_that(test$p.value >= 0.05, is_true())
})#pass



# Case 3: lower bound is correct, but upper bound is incorret
#Error: Test failed: 'ars correctly samples from truncated normal'
#Not expected: missing value where TRUE/FALSE needed
test_that("ars correctly samples from truncated normal", {
  library(truncnorm)
  x_real<-rtruncnorm(n,-10,10)
  x_ars<-ars(dnorm,n,-10,Inf)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  #compare_densities(x_real,x_ars) 
  expect_that(test$p.value >= 0.05, is_true())
})

# Case 4: WHY? different error messages with case 3.
# case 3 is -10 to Inf and case 4 is -Inf to 10
#Error: Test failed: 'ars correctly samples from truncated normal'
#Not expected: test$p.value >= 0.05 isn't true.
test_that("ars correctly samples from truncated normal", {
  library(truncnorm)
  x_real<-rtruncnorm(n,-10,10)
  x_ars<-ars(dnorm,n,-Inf,10)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  #compare_densities(x_real,x_ars) 
  expect_that(test$p.value >= 0.05, is_true())
})

# Case 5: smaller bound
test_that("ars correctly samples from truncated normal", {
  library(truncnorm)
  x_real<-rtruncnorm(n,-10,10)
  x_ars<-ars(dnorm,n,-5,5)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  #compare_densities(x_real,x_ars) 
  expect_that(test$p.value >= 0.05, is_true())
})#pass

test_that("ars correctly samples from truncated normal", {
  library(truncnorm)
  x_real<-rtruncnorm(n,-10,10)
  x_ars<-ars(dnorm,n,-1,1)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  #compare_densities(x_real,x_ars) 
  expect_that(test$p.value >= 0.05, is_true())
})#fail

################################################################################
##### Exponential dist
# Case 1: Correct lower and upper bounds
#How to find a good bound for exponential

# Case 2: No bounds
#Not expected: missing value where TRUE/FALSE needed
test_that("ars correctly samples from exponential", {
  x_real<-rexp(n)
  x_ars<-ars(dexp,n, 0)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  expect_that(test$p.value >= 0.05, is_true())
})
# Case 3: Correct upper bound but wrong lower bound
#Error: Test failed: 'ars correctly samples from exponential'
#Not expected: test$p.value >= 0.05 isn't true.
test_that("ars correctly samples from exponential", {
  x_real<-rexp(n)
  x_ars<-ars(dexp,n,1,1000)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  expect_that(test$p.value >= 0.05, is_true())
})
#Error: Test failed: 'ars correctly samples from exponential'
# "Truncated distribution: the leftmost point is the mode."
test_that("ars correctly samples from exponential", {
  x_real<-rexp(n)
  x_ars<-ars(dexp,n,0.1,1000)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  expect_that(test$p.value >= 0.05, is_true())
})
#Error: Test failed: 'ars correctly samples from exponential'
#Not expected: Bad density: not log-concave
test_that("ars correctly samples from exponential", {
  x_real<-rexp(n)
  x_ars<-ars(dexp,n,0.01,1000)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  expect_that(test$p.value >= 0.05, is_true())
})






################################################################################
##### Uniform dist
#How to solve this?
#"Uniform distribution: runif is used to generate sample"
# Case 1: Correct lower and upper bounds
test_that("ars correctly samples from uniform", {
  x_real<-runif(n,-5,5)
  x_ars<-ars(function(x) 0.1,n,-5,5)
  test <- ks.test(x_ars, x_real
  expect_that(test$p.value >= 0.05, is_true())
})
# Case 2: No bounds
test_that("ars correctly samples from uniform", {
  x_real<-runif(n,-5,5)
  x_ars<-ars(function(x) 0.1,n)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})
# Case 3: Correct lower bound but wrong upper bound
test_that("ars correctly samples from uniform", {
  x_real<-runif(n,-5,5)
  x_ars<-ars(function(x) 0.1,n,-5,Inf)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})
# Case 4: Correct upper bound but wrong lower bound
test_that("ars correctly samples from uniform", {
  x_real<-runif(n,-5,5)
  x_ars<-ars(function(x) 0.1,n,-Inf,5)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})



################################################################################
###### Laplace distribution (double exponential)
# Case 1: Correct lower and upper bounds
test_that("ars correctly samples from laplace", {
  library(smoothmest)
  laplace_pdf <- function(x) ddoublex(x, mu=0, lambda=1)
  x_real <- rdoublex(n,mu=0,lambda=1)
  x_ars <- ars(laplace_pdf,n,-10,10)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})
# Case 2: No bounds
test_that("ars correctly samples from laplace", {
  library(smoothmest)
  laplace_pdf <- function(x) ddoublex(x, mu=0, lambda=1)
  x_real <- rdoublex(n,mu=0,lambda=1)
  x_ars <- ars(laplace_pdf,n)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})#pass
# Case 3: incorrect bound
#Error: Test failed: 'ars correctly samples from laplace'
#Not expected: Bad density: not log-concave
test_that("ars correctly samples from laplace", {
  library(smoothmest)
  laplace_pdf <- function(x) ddoublex(x, mu=0, lambda=1)
  x_real <- rdoublex(n,mu=0,lambda=1)
  x_ars <- ars(laplace_pdf,n,1,2)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})
#Error: Test failed: 'ars correctly samples from laplace'
#Not expected: test$p.value >= 0.05 isn't true.
test_that("ars correctly samples from laplace", {
  library(smoothmest)
  laplace_pdf <- function(x) ddoublex(x, mu=0, lambda=1)
  x_real <- rdoublex(n,mu=0,lambda=1)
  x_ars <- ars(laplace_pdf,n,1)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})
#Not expected: test$p.value >= 0.05 isn't true.
test_that("ars correctly samples from laplace", {
  library(smoothmest)
  laplace_pdf <- function(x) ddoublex(x, mu=0, lambda=1)
  x_real <- rdoublex(n,mu=0,lambda=1)
  x_ars <- ars(laplace_pdf,n,-Inf,1)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})
#Not expected: test$p.value >= 0.05 isn't true.
test_that("ars correctly samples from laplace", {
  library(smoothmest)
  laplace_pdf <- function(x) ddoublex(x, mu=0, lambda=1)
  x_real <- rdoublex(n,mu=0,lambda=1)
  x_ars <- ars(laplace_pdf,n,0.1,0.2)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})


################################################################################
#### Test Weibull density
# Case 1: haven't found a correct bound

#Not expected: Bad density: not log-concave
test_that("ars correctly samples from weibull", {
  weibull_pdf <- function(x) dweibull(x,shape=1)
  x_real <- rweibull(n,shape=1)
  x_ars <- ars(weibull_pdf,n,0,1)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})


# Case 2: No bounds
#Not expected: missing value where TRUE/FALSE needed
test_that("ars correctly samples from weibull", {
  weibull_pdf <- function(x) dweibull(x,shape=1)
  x_real <- rweibull(n,shape=1)
  x_ars <- ars(weibull_pdf,n)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})
# Case 3: Correct lower bound but wrong upper bound
#[1] "Truncated distribution: the leftmost point is the mode."
#Error: Test failed: 'ars correctly samples from weibull'
#Not expected: test$p.value >= 0.05 isn't true.
test_that("ars correctly samples from weibull", {
  weibull_pdf <- function(x) dweibull(x,shape=1)
  x_real <- rweibull(n,shape=1)
  x_ars <- ars(weibull_pdf,n,0.0000001,1)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})
#Not expected: missing value where TRUE/FALSE needed
#Not expected: Bad density: not log-concave
test_that("ars correctly samples from weibull", {
  weibull_pdf <- function(x) dweibull(x,shape=1)
  x_real <- rweibull(n,shape=1)
  x_ars <- ars(weibull_pdf,n,0.1,1)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})
# Case 4: Correct upper bound but wrong lower bound



################################################################################
#### Test chi-square density
# Case 1: Correct lower and upper bounds
test_that("ars correctly samples from chi-square", {
  chisq_pdf <- function(x) dchisq(x,3)
  x_real <- rchisq(n,3)
  x_ars <- ars(chisq_pdf,n,0.0001,Inf) 
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})#pass
# Case 2: No bounds
#Error: Test failed: 'ars correctly samples from chi-square'
#Not expected: Density not log-concave: no maxima
test_that("ars correctly samples from chi-square", {
  chisq_pdf <- function(x) dchisq(x,3)
  x_real <- rchisq(n,3)
  x_ars <- ars(chisq_pdf,n) 
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})#pass
# Case 3: Correct lower bound but wrong upper bound

# Case 4: Correct upper bound but wrong lower bound
#Not expected: missing value where TRUE/FALSE needed
test_that("ars correctly samples from chi-square", {
  chisq_pdf <- function(x) dchisq(x,3)
  x_real <- rchisq(n,3)
  x_ars <- ars(chisq_pdf,n,0,100) 
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})


################################################################################
#### Test logistic distribution  #pass all
# Case 1: Correct lower and upper bounds
test_that("ars correctly samples from logistic", {
  x_real <- rlogis(n)
  x_ars <- ars(dlogis, n, -Inf, Inf) #-Inf, Inf
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})#pass

# Case 2: No bounds
test_that("ars correctly samples from logistic", {
  x_real <- rlogis(n)
  x_ars <- ars(dlogis, n) #-Inf, Inf
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})
# Case 3: Correct lower bound but wrong upper bound

# Case 4: Correct upper bound but wrong lower bound



################################################################################
#### Test extreme value distribution
# Case 1: Correct lower and upper bounds
test_that("ars correctly samples from extreme value", {
  library(evd)
  gev_pdf <- function(x) {dgev(x, loc = 0, scale = 1, shape = 0)}
  x_ars <- ars(gev_pdf, n, -5, 5) 
  x_real <- rgev(n, loc = 0, scale = 1, shape = 0)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})#pass
# Case 2: No bounds
#Not expected: missing value where TRUE/FALSE needed
test_that("ars correctly samples from extreme value", {
  library(evd)
  gev_pdf <- function(x) {dgev(x, loc = 0, scale = 1, shape = 0)}
  x_ars <- ars(gev_pdf, n) 
  x_real <- rgev(n, loc = 0, scale = 1, shape = 0)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})
# Case 3: Correct lower bound but wrong upper bound
#[1] "Truncated distribution: the leftmost point is the mode."
#Error: Test failed: 'ars correctly samples from extreme value'
#Not expected: test$p.value >= 0.05 isn't true.
test_that("ars correctly samples from extreme value", {
  library(evd)
  gev_pdf <- function(x) {dgev(x, loc = 0, scale = 1, shape = 0)}
  x_ars <- ars(gev_pdf, n,1) 
  x_real <- rgev(n, loc = 0, scale = 1, shape = 0)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})
# Case 4: Correct upper bound but wrong lower bound
#Not expected: missing value where TRUE/FALSE needed
test_that("ars correctly samples from extreme value", {
  library(evd)
  gev_pdf <- function(x) {dgev(x, loc = 0, scale = 1, shape = 0)}
  x_ars <- ars(gev_pdf, n,-Inf,1) 
  x_real <- rgev(n, loc = 0, scale = 1, shape = 0)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})


################################################################################
#### Test gamma distribution if the shape parameter is >= 1
# Case 1: Correct lower and upper bounds

# Case 2: No bounds

# Case 3: Correct lower bound but wrong upper bound

# Case 4: Correct upper bound but wrong lower bound

############################### Informal test ##################################
gamma_pdf <- function(x) {dgamma(x, 2)}
x_ars <- ars(gamma_pdf, n, 0.00001, Inf)
x_real <- rgamma(n, 2)
compare_densities(x_real,x_ars) #OK
ks.test(x_ars, x_real) #p-value = 0.2106


################################################################################
#### Test beta distribution if both shape parameters are >= 1
# Case 1: Correct lower and upper bounds

# Case 2: No bounds

# Case 3: Correct lower bound but wrong upper bound

# Case 4: Correct upper bound but wrong lower bound

############################### Informal test ##################################
beta_pdf <- function(x) {dbeta(x, 2, 2)}
x_ars <- ars(beta_pdf, n, 0.0001, 0.9999)
x_real <- rbeta(n, 2, 2)
compare_densities(x_real,x_ars) #OK
ks.test(x_ars, x_real) # p-value = 0.4676




################################################################################
############################### Not log-concave ################################
# Test Student's t-distribution
t_pdf <- function(x) {dt(x, df = 50)}
ars(t_pdf, n, -50, 50) # OK

# Test Cauchy distribution
cauchy_pdf <- function(x) {dcauchy(x)}
ars(cauchy_pdf, n, -Inf, Inf) # OK

# Test Pareto distribution
library(PtProcess)
pareto_pdf <- function(x) {dpareto(x, lambda = 3, a = 1)}
ars(pareto_pdf, n, 3, Inf) 

# Test F distribution
f_pdf <- function(x) {df(x, df1 = 10, df2 = 15)}
ars(f_pdf, n, 0.00001)  

# Test Chi-square with df = 1
chisq_pdf <- function(x) dchisq(x,1)
ars(chisq_pdf,n,0.001,Inf) # ok

tmp <- tempfile()
Rprof(tmp, interval = 0.01)
a <- ars(dnorm,n,-10,10)
Rprof(NULL)
summaryRprof(tmp)


# Maybe we can use the trapezoidal rule to integrate instead?
# trapezoid <- function(fun, a, b, n=100) {
#   # numerical integral of fun from a to b
#   # using the trapezoid rule with n subdivisions
#   # assume a < b and n is a positive integer
#   h <- (b-a)/n
#   x <- seq(a, b, by=h)
#   y <- fun(x)
#   s <- h * (y[1]/2 + sum(y[2:n]) + y[n+1]/2)
#   return(s)
# }
