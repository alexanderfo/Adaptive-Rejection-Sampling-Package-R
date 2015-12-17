# @knitr ars, include=FALSE
source("../ars/R/ars.R")
n <- 1000

# @knitr compare_density, include=FALSE
compare_densities<-function(x_real, x_ars, density_name) {
  p_val <- ks.test(x_ars, x_real)$p.value
  par(mfrow=c(1,2))
  plot(density(x_real),type="l", 
       main=paste(density_name, " v.s. ARS:\n p-value ", round(p_val,4), sep = ""),
       xlab="x", lty = 2, lwd = 2)
  lines(density(x_ars),type="l", lty = 1, lwd = 2)
  legend(x = "topright", legend = c(paste("True ", density_name, sep=""), "ARS"), 
         cex = 0.8, lty = c(2, 1), lwd = c(2, 2))
  qqplot(x_ars, x_real, main=paste("QQplot\n ", density_name, " v.s. ARS", sep = ""),
         xlab = "True", ylab = "ARS samples", pch = 16)
  par(mfrow=c(1,1))
}

# @knitr normal, echo = FALSE
set.seed(112)
x_ars<-ars(dnorm,1000)
x_real<-rnorm(1000)
compare_densities(x_real, x_ars, "Normal")

# @knitr truncated_normal, echo = FALSE
set.seed(111)
x_ars<-ars(dnorm,1000,-10,2)
x_real<-rtruncnorm(1000,-10,2)
compare_densities(x_real, x_ars, "Truncated Normal(0, 1) from -10 to 2")

# @knitr uniform, echo = FALSE
set.seed(112)
x_real<-runif(n,-5,5)
x_ars<-ars(function(x) 0.1*x^0,n,-5,5)
compare_densities(x_real, x_ars, "Uniform(-5, 5)")

# @knitr exponential, echo = FALSE
set.seed(112)
x_real<-rexp(n)
x_ars<-ars(dexp,n, 0)
compare_densities(x_real, x_ars, "Exponential(1)")

# @knitr gamma, echo = FALSE
set.seed(7)
gamma_pdf <- function(x) {dgamma(x, 2)}
x_ars <- ars(gamma_pdf, n, 0.00001, Inf)
x_real <- rgamma(n, 2)
compare_densities(x_real, x_ars, "Gamma(2, 1)")

# @knitr beta, echo = FALSE
set.seed(7)
beta_pdf <- function(x) {dbeta(x, 2, 2)}
x_ars <- ars(beta_pdf, n, 0.0001, 0.9999)
x_real <- rbeta(n, 2, 2)
compare_densities(x_real, x_ars, "Beta(2, 2)")

# @knitr logistic, echo = FALSE
set.seed(5)
x_ars <- ars(dlogis, n, -Inf, Inf)
x_real <- rlogis(n)
compare_densities(x_real, x_ars, "Logistic")

# @knitr extreme_value, echo = FALSE, messages = FALSE
set.seed(16)
gev_pdf <- function(x) {dgev(x, loc = 0, scale = 1, shape = 0)}
x_ars <- ars(gev_pdf, n, -5, 5) #-Inf, Inf
x_real <- rgev(n, loc = 0, scale = 1, shape = 0)
compare_densities(x_real, x_ars, "EVD")

# @knitr laplace, echo = FALSE, messages = FALSE
set.seed(113)
laplace_pdf <- function(x) ddoublex(x, mu=0, lambda=1)
x_real <- rdoublex(n,mu=0,lambda=1)
x_ars <- ars(laplace_pdf,n)
compare_densities(x_real, x_ars, "Laplace")

# @knitr chisq_3, echo = FALSE
set.seed(2)
chisq_pdf <- function(x) dchisq(x,3)
x_real <- rchisq(n,3)
x_ars <- ars(chisq_pdf,n,0.0001,Inf)
compare_densities(x_real, x_ars, "Chi-squared 3")

# @knitr weibull, echo = FALSE
set.seed(115)
weibull_pdf <- function(x) dweibull(x,shape=1)
x_real <- rweibull(n,shape=1)
x_ars <- ars(weibull_pdf,n,0,Inf)
compare_densities(x_real, x_ars, "Weibull")



# @knitr chisq_1, echo = FALSE
set.seed(115)
chisq_pdf <- function(x) dchisq(x,1)
ars(chisq_pdf,n,0.001,Inf)

# @knitr student_t, echo = FALSE
set.seed(115)
t_pdf <- function(x) {dt(x, df = 50)}
ars(t_pdf, n, -50, 50)

# @knitr cauchy, echo = FALSE
set.seed(115)
cauchy_pdf <- function(x) {dcauchy(x)}
ars(cauchy_pdf, n, -Inf, Inf)

# @knitr pareto, echo = FALSE
set.seed(115)
pareto_pdf <- function(x) {dpareto(x, lambda = 3, a = 1)}
ars(pareto_pdf, n, 3, Inf)

# @knitr fdist, echo = FALSE
set.seed(115)
f_pdf <- function(x) {df(x, df1 = 10, df2 = 15)}
ars(f_pdf, n, 0.00001)

# @knitr diverging_density, echo = FALSE
dvg_pdf <- function(x) {exp(x)}
ars(dvg_pdf, n, 1, Inf)
