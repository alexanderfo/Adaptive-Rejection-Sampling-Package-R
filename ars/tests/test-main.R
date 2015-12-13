# Test the arSampler
# setwd("~/git/stat243-project/")
source("ars/R/new_ars.R")
library(testthat)

compare_densities<-function(x_real,x_ars) {
  print(density(x_real))
  print(density(x_ars))
  
  par(mfrow=c(1,1))
  plot(density(x_real),type="l",col="green",main="real dens green, approx dens red",xlab="x")
  lines(density(x_ars),type="l",col="red")
}

n<-100 #number of samples
################################################################################
############################ Log Concave Dist ##################################
#### Normal
# Case 1: Correct lower and upper bounds
# To Huohuo: this is a reference code that might work
test_that("ars correctly samples from truncated normal", {
  library(truncnorm)
  x_real<-rtruncnorm(n,a=-10,b=10)
  x_ars<-arSampler(dnorm,n,-10,10)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  expect_that(test$p.value >= 0.05, is_true())
})

# Case 2: No bounds
test_that("ars correctly samples from standard normal", {
  x_real<-rnorm(n)
  x_ars<-arSampler(dnorm,n)
  test <- ks.test(x_ars, x_real) # p-value = 0.8127
  compare_densities(x_real,x_ars) 
  expect_that(test$p.value >= 0.05, is_true())
})


# Case 3: Correct lower bound but wrong upper bound

# Case 4: Correct upper bound but wrong lower bound

############################### Informal test ##################################
x_real <- rnorm(n, 0, 1)
x_ars <- arSampler(dnorm, n, -Inf, Inf)
compare_densities(x_real,x_ars) #OK
ks.test(x_ars, x_real)
shapiro.test(x_ars)


################################################################################
##### Truncated normal

############################### Informal test ##################################

# library(truncnorm)
# x_real<-rtruncnorm(n,a=-10,b=10)
# x_ars<-arSampler(dnorm,n,-10,10)
# 
# compare_densities(x_real,x_ars) #OK
# test <- ks.test(x_ars, x_real) # p-value = 0.8127
# expect_that(test$p.value >= 0.05, is_true())


################################################################################
##### Exponential dist
# Case 1: Correct lower and upper bounds

# Case 2: No bounds

# Case 3: Correct lower bound but wrong upper bound

# Case 4: Correct upper bound but wrong lower bound


############################### Informal test ##################################
# problem with the limits can't choose upper>800, lower<0 (where the distribution 
# is 0 only numerically)
x_real <- rexp(n)
x_ars<-arSampler(dexp,n,0,Inf) 
compare_densities(x_real,x_ars) #OK
ks.test(x_ars, x_real)  # OK, p-value = 0.6994

test_that("ars correctly samples from exponential", {
  x_real <- rexp(n)
  x_ars<-arSampler(dexp,n,0,700) 
  compare_densities(x_real,x_ars) #OK
  test <- ks.test(x_ars, x_real) # p-value = 0.6994
  expect_that(test$p.value >= 0.05, is_true())
})


################################################################################
##### Uniform dist
# Case 1: Correct lower and upper bounds

# Case 2: No bounds

# Case 3: Correct lower bound but wrong upper bound

# Case 4: Correct upper bound but wrong lower bound

############################### Informal test ##################################
x_real <- runif(n,-5,5)
unif_pdf <- function(x) {dunif(x, -5, 5)}
x_ars <- arSampler(unif_pdf,n,-5,5)
compare_densities(x_real,x_ars) #OK
ks.test(x_ars, x_real) # p-value = 0.281


################################################################################
###### Laplace distribution (double exponential)
# Case 1: Correct lower and upper bounds

# Case 2: No bounds

# Case 3: Correct lower bound but wrong upper bound

# Case 4: Correct upper bound but wrong lower bound

############################### Informal test ##################################
library(smoothmest)
laplace_pdf <- function(x) ddoublex(x, mu=0, lambda=1)
x_real <- rdoublex(n,mu=0,lambda=1)
x_ars <- arSampler(laplace_pdf,n,-10,10)
compare_densities(x_real,x_ars) #OK
ks.test(x_ars, x_real) # p-value = 0.8127


################################################################################
#### Test Weibull density
# Case 1: Correct lower and upper bounds

# Case 2: No bounds

# Case 3: Correct lower bound but wrong upper bound

# Case 4: Correct upper bound but wrong lower bound

############################### Informal test ##################################
weibull_pdf <- function(x) dweibull(x,shape=1)
x_real <- rweibull(n,shape=1)
x_ars <- arSampler(weibull_pdf,n,0,Inf) #does not work for lower=0 and upper>800
compare_densities(x_real,x_ars) #OK
ks.test(x_ars, x_real) # p-value = 0.4676


################################################################################
#### Test chi-square density
# Case 1: Correct lower and upper bounds

# Case 2: No bounds

# Case 3: Correct lower bound but wrong upper bound

# Case 4: Correct upper bound but wrong lower bound

############################### Informal test ##################################
#Ok logconcave when df>=2
chisq_pdf <- function(x) dchisq(x,3)

x_real <- rchisq(n,3)
x_ars <- arSampler(chisq_pdf,n,0.0001,Inf)  # check when lb = 0

compare_densities(x_real,x_ars) #OK
ks.test(x_ars, x_real) # p-value = 0.3667


################################################################################
#### Test logistic distribution
# Case 1: Correct lower and upper bounds

# Case 2: No bounds

# Case 3: Correct lower bound but wrong upper bound

# Case 4: Correct upper bound but wrong lower bound

############################### Informal test ##################################
x_ars <- arSampler(dlogis, n, -Inf, Inf) #-Inf, Inf
x_real <- rlogis(n)
compare_densities(x_real,x_ars) #OK
ks.test(x_ars, x_real) #p-value = 0.5806


################################################################################
#### Test extreme value distribution
# Case 1: Correct lower and upper bounds

# Case 2: No bounds

# Case 3: Correct lower bound but wrong upper bound

# Case 4: Correct upper bound but wrong lower bound

############################### Informal test ##################################
library(evd)
gev_pdf <- function(x) {dgev(x, loc = 0, scale = 1, shape = 0)}
x_ars <- arSampler(gev_pdf, n, -5, 5) #-Inf, Inf
x_real <- rgev(n, loc = 0, scale = 1, shape = 0)
compare_densities(x_real,x_ars) # OK
ks.test(x_ars, x_real) #p-value = 0.9062


################################################################################
#### Test gamma distribution if the shape parameter is >= 1
# Case 1: Correct lower and upper bounds

# Case 2: No bounds

# Case 3: Correct lower bound but wrong upper bound

# Case 4: Correct upper bound but wrong lower bound

############################### Informal test ##################################
gamma_pdf <- function(x) {dgamma(x, 2)}
x_ars <- arSampler(gamma_pdf, n, 0.00001, Inf)
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
x_ars <- arSampler(beta_pdf, n, 0.0001, 0.9999)
x_real <- rbeta(n, 2, 2)
compare_densities(x_real,x_ars) #OK
ks.test(x_ars, x_real) # p-value = 0.4676




################################################################################
############################### Not log-concave ################################
# Test Student's t-distribution
t_pdf <- function(x) {dt(x, df = 50)}
arSampler(t_pdf, n, -50, 50) # OK

# Test Cauchy distribution
cauchy_pdf <- function(x) {dcauchy(x)}
arSampler(cauchy_pdf, n, -100, 100) # OK

# Test Pareto distribution
library(PtProcess)
pareto_pdf <- function(x) {dpareto(x, lambda = 3, a = 2)}
arSampler(pareto_pdf, n, 3, 100) # OK

# Test F distribution
f_pdf <- function(x) {df(x, df1 = 10, df2 = 10)}
arSampler(f_pdf, n, 0.00001, 10)  # OK

# Test Chi-square with df = 1
chisq_pdf <- function(x) dchisq(x,1)
arSampler(chisq_pdf,n,0.001,100)


tmp <- tempfile()
Rprof(tmp, interval = 0.01)
a <- arSampler(dnorm,n,-10,10)
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
