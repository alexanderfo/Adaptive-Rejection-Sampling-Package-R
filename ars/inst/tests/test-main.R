library(testthat)
library(truncnorm) # for truncated normal
library(smoothmest) # for laplace
library(PtProcess) # for pareto
library(truncdist) # for all kinds of truncations

cat(paste (' RUNNING TESTS! n=1000 samples (dot indicates passed)', '\n'))

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
  a<-expect_that(test$p.value >= 0.01, is_true())
  print("Test: Normal dist, bounds=(-Inf,Inf)")
})

# Case 2: Correct finite bounds
test_that("ars correctly samples from truncated normal", {
  
  library(truncnorm)
  x_real<-rtruncnorm(n,a=-10,b=10)
  x_ars<-ars(dnorm,n,-10,10)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Normal dist, truncated, bounds=[-10,10]")
})

# Case 3: Infinite lb, finite ub
test_that("ars correctly samples from truncated normal where lb is infinite", {
  
  library(truncnorm)
  x_real<-rtruncnorm(n,a=-Inf,b=10)
  x_ars<-ars(dnorm, n, -Inf, 10)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Normal dist, truncated, bounds=[-Inf,10]")
})

# Case 4: finite lb, infinite ub
test_that("ars correctly samples from truncated normal where ub is infinite", {
  
  library(truncnorm)
  x_real<-rtruncnorm(n,a=2,b=Inf)
  x_ars<-ars(dnorm, n, 2, Inf)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Normal dist, truncated, bounds=[2,Inf]")
})

# Case 5: scaled dnorm
test_that("ars correctly samples from scaled standard normal", {
  f <- function(x) 0.5*dnorm(x)
  x_real<-rnorm(n)
  x_ars<-ars(f,n)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Normal dist scaled, 0.5*dnorm, bounds=(-Inf,Inf)")
})

# Case 6: shifted dnorm
test_that("ars correctly samples from shifted standard normal", {
  
  f <- function(x) dnorm(x,mean=5)
  x_real<-rnorm(n)+5
  x_ars<-ars(f,n,-100,100)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Normal dist shifted, mean=5, bounds=(-Inf,Inf)")
})

################################################################################
# 2. Exponential dist

# Case 1: No bounds
test_that("ars correctly samples from exponential with infinite ub", {
  
  x_real<-rexp(n)
  x_ars<-ars(dexp,n, 0)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Exponential dist, bounds=[0,Inf]")
})

# Case 2: Finite ub
test_that("ars correctly samples from exponential with finite ub", {
  
  x_real<-rexp(n)
  x_ars<-ars(dexp, n, 0, 20)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Exponential dist, truncated, bounds=[0,20]")
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
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Uniform dist, bounds=[0,1]")
})

# Case 2: Correct lower and upper bounds other than [0,1]
test_that("ars correctly samples from uniform(-5,5)", {
  
  x_real<-runif(n,-5,5)
  x_ars<-ars(function(x) 0.1*x^0, n, -5, 5)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Uniform dist, bounds=[-5,5]")
})

test_that("ars correctly samples from uniform (13, 27)", {
  
  x_real<-runif(n, 13, 27)
  x_ars<-ars(function(x) 0.1*x^0, n, 13, 27)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Uniform dist, bounds=[13,27]")
})

# Case 3: Infinite boundaries - should diverge
test_that("ars correctly report error when infinite bound is given to uniform", {
  expect_error(x_ars<-ars(function(x) 0.1,n), "Bad density: density diverges")
  expect_error(x_ars<-ars(function(x) 0.1, n, -5, Inf), "Bad density: density diverges")
  expect_error(x_ars<-ars(function(x) 0.1, n, -Inf, 37), "Bad density: density diverges")
  expect_error(x_ars<-ars(function(x) 17,n), "Bad density: density diverges")
})
print("4 Tests: Four different constant functions on bounds=(-Inf,Inf) (should diverge)")

################################################################################
# 4. Test gamma distribution if the shape parameter is >= 1
# Case 1: Correct lower and upper bounds
test_that("ars correctly sample from gamma(2,1)", {
  
  gamma_pdf <- function(x) {dgamma(x, 2)}
  x_real <- rgamma(n, 2)
  x_ars <- ars(gamma_pdf, n, 0.00001, Inf)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.01, is_true())
})
print("Test: Gamma dist, shape=2, bounds=(-Inf,Inf)")

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
  x_ars <- ars(beta_pdf, n, 0.000001, 0.9999999)
  x_real <- rbeta(n, shape1, shape2)
  test <- ks.test(x_ars, x_real) # p-value = 0.4676
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Beta dist, shape1=1, shape2=1, bounds=(-Inf,Inf)")

  shape1 <- 2
  shape2 <- 1
  beta_pdf <- function(x) {dbeta(x, shape1, shape2)}
  x_ars <- ars(beta_pdf, n, 0.000001, 0.9999999)
  x_real <- rbeta(n, shape1, shape2)
  test <- ks.test(x_ars, x_real) # p-value = 0.4676
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Beta dist, shape1=2, shape2=1, bounds=(-Inf,Inf)")

  shape1 <- 2
  shape2 <- 3
  beta_pdf <- function(x) {dbeta(x, shape1, shape2)}
  x_ars <- ars(beta_pdf, n, 0.000001, 0.9999999)
  x_real <- rbeta(n, shape1, shape2)
  test <- ks.test(x_ars, x_real) # p-value = 0.4676
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Beta dist, shape1=2, shape2=3, bounds=(-Inf,Inf)")
})


################################################################################
# 6. Test logistic distribution  #pass all
# Case 1: No bound
test_that("ars correctly samples from logistic (0,1)", {
  x_ars <- ars(dlogis, n, -Inf, Inf) #-Inf, Inf
  x_real <- rlogis(n)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Logistic dist, location=0, scale=1, bounds=(-Inf,Inf)")
})#pass

test_that("ars correctly samples from logistic (2,2)", {
  set.seed(0) #unstable in very rare cases
  location <- 2
  scale <- 2
  f <- function(x) dlogis(x, location, scale)
  x_ars <- ars(f, n, -Inf, Inf) #-Inf, Inf
  x_real <- rlogis(n, location, scale)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Logistic dist, location=2, scale=2, bounds=(-Inf,Inf)")
})#pass

################################################################################
# 7. Test extreme value distribution
# Case 1: Truncated range
test_that("ars correctly samples from extreme value dist with truncated range", {
  
  gev_pdf <- function(x) {dgev(x, loc = 0, scale = 1, shape = 0)}
  x_ars <- ars(gev_pdf, n, -5, 5)
  x_real <- rtrunc(n, "gev", a=-5, b=5)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Extreme Value dist, truncated, location=0, scale=1, shape=0, bounds=[-5,5]")
})#pass

# Case 2: No bounds
test_that("ars correctly samples from extreme value", {
  
  gev_pdf <- function(x) {dgev(x, loc = 0, scale = 1, shape = 0)}
  x_ars <- ars(gev_pdf, n,-500,500)
  x_real <- rgev(n, loc = 0, scale = 1, shape = 0)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Extreme Value dist, location=0, scale=1, shape=0, bounds=(-Inf,Inf)")
})

################################################################################
# 8. Laplace distribution (double exponential)

# Case 1: No bounds
test_that("ars correctly samples from laplace with infinite range", {
  
  laplace_pdf <- function(x) ddoublex(x, mu=0, lambda=1)
  x_real <- rdoublex(n,mu=0,lambda=1)
  x_ars <- ars(laplace_pdf,n)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Laplace dist, location=0, scale=1, shape=0, bounds=(-Inf,Inf)")
})#pass

################################################################################
# 9. Test chi-square density
# Case 1: Infinite ub
test_that("ars correctly samples from chi-square with infinite ub", {
  chisq_pdf <- function(x) dchisq(x,3)
  x_real <- rchisq(n,3)
  x_ars <- ars(chisq_pdf,n,0.0001,Inf)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Chi-Squared dist, df=3, bounds=[0,Inf)")
})

# Case 2: Finite ub
test_that("ars correctly samples from chi-square with finite ub", {
  
  chisq_pdf <- function(x) dchisq(x,3)
  x_real <- rtrunc(n, "chisq", a=0, b=50,df=3)
  x_ars <- ars(chisq_pdf,n,0.01,50)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Chi-Squared dist, truncated, df=3, bounds=[0,50]")
})#pass

# Case 3: Truncated
test_that("ars correctly samples from chi-square with truncated range", {
  chisq_pdf <- function(x) dchisq(x,3)
  x_real <- rtrunc(n, "chisq", a=12, b=50,df=3)
  x_ars <- ars(chisq_pdf,n,12,50)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Chi-Squared dist, truncated, df=3, bounds=[12,50]")
})#pass


################################################################################
# 10. Test Weibull density
# Case 1: haven't found a correct bound

test_that("ars correctly samples from weibull with finite ub", {
  
  weibull_pdf <- function(x) dweibull(x,shape=1)
  x_real <- rtrunc(n, "weibull", a=0, b=5, shape=1)
  x_ars <- ars(weibull_pdf,n,0,5)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Weibull dist, truncated, shape=1, bounds=[0,5]")
})

# Case 2: Infinite ub
test_that("ars correctly samples from weibull with infinite ub", {
  
  weibull_pdf <- function(x) dweibull(x,shape=1)
  x_real <- rweibull(n,shape=1)
  x_ars <- ars(weibull_pdf,n,0)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Weibull dist, shape=1, bounds=[0,Inf)")
})

# Case 3: Truncated range
test_that("ars correctly samples from weibull with truncated range", {
  
  weibull_pdf <- function(x) dweibull(x,shape=1)
  x_real <- rtrunc(n, "weibull", a=5, b=10, shape=1)
  x_ars <- ars(weibull_pdf,n,5,10)
  test <- ks.test(x_ars, x_real)
  expect_that(test$p.value >= 0.01, is_true())
  print("Test: Weibull dist truncated, shape=1, bounds=[5,10]")
})


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
  print("Test that Chi-sq dist, df=1 is not log-concave")
  expect_error(ars(t_pdf, n, -50, 50), "not log-concave")
  print("Test that Student t dist, df=50 is not log-concave")
  expect_error(ars(cauchy_pdf, n, -Inf, Inf), "not log-concave")
  print("Test that Cauchy dist is not log-concave")
  expect_error(ars(pareto_pdf, n, 3, Inf), "not log-concave")
  print("Test that Pareto dist shape=3 is not log-concave")
  expect_error(ars(f_pdf, n, 0.00001), "not log-concave")
  print("Test that F dist df1=10, df2=15 is not log-concave")
})
cat(paste ('\n', '35 tests completed! (dot indicates test passed)','\n'))


