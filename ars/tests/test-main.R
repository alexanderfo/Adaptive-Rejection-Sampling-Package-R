# Test the arSampler
setwd("~/git/stat243-project/")
source("ars/R/ars.R")
library(truncnorm)

compare_densities<-function(x_real,x_ars) {
  print(density(x_real))
  print(density(x_ars))
  
  par(mfrow=c(1,1))
  plot(density(x_real),type="l",col="green",main="real dens green, approx dens red",xlab="x")
  lines(density(x_ars),type="l",col="red")
}


n<-500 #number of samples


##### Truncated normal
x_real<-rtruncnorm(n,a=-10,b=10)
x_ars<-arSampler(dnorm,n,-10,10)

compare_densities(x_real,x_ars) #OK



##### Exponential dist, PARTIALLY DOES NOT WORK
x_real <- rexp(n)
x_ars<-arSampler(dexp,n,0.000001,700) #problem with the limits can't choose upper>800, lower=0
#Problem in draw_sample with the integrate function

compare_densities(x_real,x_ars)



##### Uniform dist, DOES NOT WORK
x_real <- runif(n)
x_ars <- arSampler(dunif,n,0.5,0.6)
#Problem in draw_sample with the integrate function


###### Laplace distribution (double exponential), DOES NOT WORK 

# in draw sample integrate does not work: numerator <- sapply(1:length(envelope), function(i) integrate(envelope[[i]], z[i], z[i+1])$value)
# Produces two errors: Error in integrate(envelope[[i]], z[i], z[i + 1]) : a limit is missing (NaN is created)
# and Error in integrate(envelope[[i]], z[i], z[i + 1]) : non-finite function value
library(smoothmest)

laplace_pdf <- function(x) ddoublex(x, mu=0, lambda=1)

x_real <- rdoublex(n,mu=0,lambda=1)
x_ars <- arSampler(laplace_pdf,n,-1,1) #does not work with a greater interval if we include zero in the interval
# I changed in draw_sample, then it worked for some intervals, very unstable

compare_densities(x_real,x_ars)


# Test Weibull density

weibull_pdf <- function(x) dweibull(x,shape=1)

x_real <- rweibull(n,shape=1)
x_ars <- arSampler(weibull_pdf,n,0.0001,700) #does not work for lower=0 and upper>800


compare_densities(x_real,x_ars) #OK


# Test chi-square density DOES NOT WORK

chisq_pdf <- function(x) dchisq(x,1)
x_ars <- arSampler(chisq_pdf,n,0.00001,100) # is not log-concave for parameter, degrees of freedom = 1.
is_logconcave(chisq_pdf,0.00001,100)

#Ok
chisq_pdf <- function(x) dchisq(x,2)
is_logconcave(chisq_pdf,0.00001,100)
#Ok logconcave when df>=2

x_real <- rchisq(n,2)
x_ars <- arSampler(chisq_pdf,n,0.00001,100)

compare_densities(x_real,x_ars) #OK






# Maybe we can use the trapezoidal rule to integrate instead?
trapezoid <- function(fun, a, b, n=100) {
  # numerical integral of fun from a to b
  # using the trapezoid rule with n subdivisions
  # assume a < b and n is a positive integer
  h <- (b-a)/n
  x <- seq(a, b, by=h)
  y <- fun(x)
  s <- h * (y[1]/2 + sum(y[2:n]) + y[n+1]/2)
  return(s)
}
