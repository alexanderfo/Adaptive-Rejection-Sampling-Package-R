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


n<-100 #number of samples


##### Truncated normal
x_real<-rtruncnorm(n,a=-10,b=10)
x_ars<-arSampler(dnorm,n,-10,10)

compare_densities(x_real,x_ars) #OK



##### Exponential dist
#problem with the limits can't choose upper>800, lower<0 (where the distribution is 0 only numerically)
x_real <- rexp(n)
x_ars<-arSampler(dexp,n,0,700) 
compare_densities(x_real,x_ars) #OK



##### Uniform dist
x_real <- runif(n,-5,5)
x_ars <- arSampler(dunif,n,-5,5, -5, 5)
compare_densities(x_real,x_ars) #OK


###### Laplace distribution (double exponential)
library(smoothmest)

laplace_pdf <- function(x) ddoublex(x, mu=0, lambda=1)

x_real <- rdoublex(n,mu=0,lambda=1)
x_ars <- arSampler(laplace_pdf,n,-10,10)

compare_densities(x_real,x_ars) #OK


# Test Weibull density

weibull_pdf <- function(x) dweibull(x,shape=1)

x_real <- rweibull(n,shape=1)
x_ars <- arSampler(weibull_pdf,n,0,700) #does not work for lower=0 and upper>800


compare_densities(x_real,x_ars) #OK


# Test chi-square density

chisq_pdf <- function(x) dchisq(x,1)
x_ars <- arSampler(chisq_pdf,n,0.00001,100) # is not log-concave for parameter, degrees of freedom = 1.


#Ok logconcave when df>=2
chisq_pdf <- function(x) dchisq(x,3)

x_real <- rchisq(n,3)
x_ars <- arSampler(chisq_pdf,n,0,100)

compare_densities(x_real,x_ars) #OK


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
