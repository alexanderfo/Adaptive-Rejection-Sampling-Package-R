# Test the ars
# setwd("~/git/stat243-project/")
source("ars/R/ars.R")
library(testthat)
library(truncnorm) # for truncated normal
library(smoothmest) # for laplace
library(PtProcess) # for pareto

# number of samples
n<-1000

################################################################################
############################ Log Concave Dist ##################################
# 1. Normal
# Case 1: No bounds
test_that("ars correctly samples from standard normal with infinte bounds", {
  x_real<-rnorm(n)
  x_ars<-ars(dnorm,n)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  expect_that(test$p.value >= 0.05, is_true())
})

# Case 2: Correct finite bounds
test_that("ars correctly samples from truncated normal", {
  library(truncnorm)
  x_real<-rtruncnorm(n,a=-10,b=10)
  x_ars<-ars(dnorm,n,-10,10)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  expect_that(test$p.value >= 0.05, is_true())
})

# Case 3: Infinite lb, finite ub
test_that("ars correctly samples from truncated normal where lb is infinite", {
  library(truncnorm)
  x_real<-rtruncnorm(n,a=-Inf,b=10)
  x_ars<-ars(dnorm, n, -Inf, 10)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  expect_that(test$p.value >= 0.05, is_true())
})

# Case 4: finite lb, infinite ub
test_that("ars correctly samples from truncated normal where ub is infinite", {
  library(truncnorm)
  x_real<-rtruncnorm(n,a=2,b=Inf)
  x_ars<-ars(dnorm, n, 2, Inf)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  expect_that(test$p.value >= 0.05, is_true())
})

# Case 5: scaled dnorm
test_that("ars correctly samples from standard normal with infinte bounds", {
  f <- function(x) 0.5*dnorm(x)
  x_real<-rnorm(n)
  x_ars<-ars(f,n)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  compare_densities(x_real, x_ars)
  expect_that(test$p.value >= 0.05, is_true())
})

################################################################################
# 2. Exponential dist

# Case 1: No bounds
test_that("ars correctly samples from exponential with infinite ub", {
  x_real<-rexp(n)
  x_ars<-ars(dexp,n, 0)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  expect_that(test$p.value >= 0.05, is_true())
})

# Case 2: Finite ub
test_that("ars correctly samples from exponential with finite ub", {
  x_real<-rexp(n)
  x_ars<-ars(dexp, n, 0, 200)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  expect_that(test$p.value >= 0.05, is_true())
})

################################################################################
# 3. Uniform dist
#How to solve this? DONE - problem with R's integrate function
#"Uniform distribution: runif is used to generate sample"

# Case 1: dunif
test_that("ars correctly samples from uniform(0,1)", {
  x_real<-runif(n)
  x_ars<-ars(dunif, n, 0, 1)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})

# Case 2: Correct lower and upper bounds other than [0,1]
test_that("ars correctly samples from uniform(-5,5)", {
  x_real<-runif(n,-5,5)
  x_ars<-ars(function(x) 0.1*x^0, n, -5, 5)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})

test_that("ars correctly samples from uniform (13, 27)", {
  x_real<-runif(n, 13, 27)
  x_ars<-ars(function(x) 0.1*x^0, n, 13, 27)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})

# Case 3: Infinite boundaries - should diverge
test_that("ars correctly report error when infinite bound is given to uniform", {
  expect_error(x_ars<-ars(function(x) 0.1,n), "Bad density: density diverges")
  expect_error(x_ars<-ars(function(x) 0.1, n, -5, Inf), "Bad density: density diverges")
  expect_error(x_ars<-ars(function(x) 0.1, n, -Inf, 37), "Bad density: density diverges")
  expect_error(x_ars<-ars(function(x) 17,n), "Bad density: density diverges")
})


################################################################################
# 4. Test gamma distribution if the shape parameter is >= 1
# Case 1: Correct lower and upper bounds
test_that("ars correctly sample from gamma(2,1)", {
  gamma_pdf <- function(x) {dgamma(x, 2)}
  x_real <- rgamma(n, 2)
  x_ars <- ars(gamma_pdf, n, 0.00001, Inf)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})

# Case 2: No bounds

# Case 3: Correct lower bound but wrong upper bound

# Case 4: Correct upper bound but wrong lower bound


################################################################################
# 5. Test beta distribution if both shape parameters are >= 1
# Case 1: Correct lower and upper bounds
test_that("ars correctly sample from beta with both shape >= 1", {
  shape1 <- 1
  shape2 <- 1
  beta_pdf <- function(x) {dbeta(x, shape1, shape2)}
  x_ars <- ars(beta_pdf, n, 0.0001, 0.99999)
  x_real <- rbeta(n, shape1, shape2)
  test <- ks.test(x_ars, x_real) # p-value = 0.4676
  expect_that(test$p.value >= 0.05, is_true())
  
  shape1 <- 2
  shape2 <- 1
  beta_pdf <- function(x) {dbeta(x, shape1, shape2)}
  x_ars <- ars(beta_pdf, n, 0.0001, 0.99999)
  x_real <- rbeta(n, shape1, shape2)
  test <- ks.test(x_ars, x_real) # p-value = 0.4676
  expect_that(test$p.value >= 0.05, is_true())
  
  shape1 <- 2
  shape2 <- 3
  beta_pdf <- function(x) {dbeta(x, shape1, shape2)}
  x_ars <- ars(beta_pdf, n, 0.0001, 0.99999)
  x_real <- rbeta(n, shape1, shape2)
  test <- ks.test(x_ars, x_real) # p-value = 0.4676
  expect_that(test$p.value >= 0.05, is_true())
})


################################################################################
# 6. Test logistic distribution  #pass all
# Case 1: No bound
test_that("ars correctly samples from logistic (0,1)", {
  set.seed(5)
  x_ars <- ars(dlogis, n, -Inf, Inf) #-Inf, Inf
  x_real <- rlogis(n)
  compare_densities(x_real, x_ars)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})#pass

test_that("ars correctly samples from logistic (2,2)", {
  set.seed(5)
  location <- 2
  scale <- 2
  f <- function(x) dlogis(x, location, scale)
  x_ars <- ars(f, n, -Inf, Inf) #-Inf, Inf
  x_real <- rlogis(n, location, scale)
  compare_densities(x_real, x_ars)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.05, is_true())
})#pass

################################################################################
# 7. Test extreme value distribution
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


################################################################################
# 8. Laplace distribution (double exponential)
# Case 1: Correct lower and upper bounds
test_that("ars correctly samples from laplace", {
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
# 9. Test chi-square density
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
# 10. Test Weibull density
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
############################### Not log-concave ################################
# Test Chi-square with df = 1
test_that("ars correctly report error on non log-concavity", {
  chisq_pdf <- function(x) dchisq(x,1)
  t_pdf <- function(x) dt(x, df = 50) 
  cauchy_pdf <- function(x) dcauchy(x)
  pareto_pdf <- function(x) dpareto(x, lambda = 3, a = 1)
  f_pdf <- function(x) df(x, df1 = 10, df2 = 15)
  
  expect_error(ars(chisq_pdf,n,0.001,Inf), "not log-concave")
  expect_error(ars(t_pdf, n, -50, 50), "not log-concave") 
  expect_error(ars(cauchy_pdf, n, -Inf, Inf), "not log-concave")
  expect_error(ars(pareto_pdf, n, 3, Inf), "not log-concave")
  expect_error(ars(f_pdf, n, 0.00001), "not log-concave")
})

