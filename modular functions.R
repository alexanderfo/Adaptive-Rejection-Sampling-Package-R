# NAMING CONVENTION!

# noun for variables
# verb for functions
# underscore to link words
# lowercase for everything

matrix of vertices and derivatives(slopes)

create_segment <- function(v1, v2){

}


# Ji and Alex
update_u <- fucntion(vertices, slopes){
	pass
}

update_l <- fucntion(vertices){
	pass
}

# Mengfei & Daisy
draw_sample <- function(envelope, z, num_of_samples = 1){
	# Draw sampling point candidates
  # Args: 
  #   envelope: vector of exponentiated linear envelope functions
  #   z: vector of intersection points
  # Output:
  #   candidates: vector sampling candidates
  numerator <- sapply(1:length(envelope), function(i) integrate(envelope[[i]], z[i], z[i+1])$value)
  denominator <- sum(numerator)
  cdf_end_pts <- cumsum(numerator/denominator)
  w <- runif(num_of_samples)
  # locate the bin each w falls in, i.e. find i
  w_idx <- sapply(w, function(i) which(i <= cdf_end_pts)[1])
  f <- function(x, i, w_i) {
    # Args:
    #   x: one sampling candidate
    #   i: the index of the Uniform variable in the w vector
    #   w_i: the bin index that the Uniform variable falls in
    integrate(envelope[[w_i]], lower = z[w_i], upper = x)$value/denominator - w[i] + ifelse(w_i == 1, 0, cdf_end_pts[w_i-1])
    }
  candidates <- sapply(1:num_of_samples, 
                       function(i) {
                         uniroot(f,interval = c(z[w_idx[i]], z[w_idx[i]+1]), 
                                 i, w_idx[i], tol = 1e-6)$root})
  return(candidates)
}

# little test
f1 <- function(x) exp(x)
f2 <- function(x) exp(-(x-1))
envelope <- c(f1, f2)
z <- c(-1, 0.5, 2)
# see code efficiency
tmp <- tempfile()
Rprof(tmp, interval = 0.1)
a=draw_sample(envelope, z, num_of_samples = 100000)
Rprof(NULL)
summaryRprof(tmp)


insert_new_vertex(vertices, new_vertex){
	pass
}

evaluate_deriv(h){
	pass
}

create_z(u){
	return(z)
}

# Mengfei, Zhuangdi
f <- function(x) {return(x^(2-1) * (1-x)^(2-1))}
# x must be between [0,1]
check_support_boundaries <- function(f, lower, upper) {
  f_lval <- f(lower)
  f_uval <- f(upper)
  # check legitimacy of density
  if (f_lval == Inf | f_lval == -Inf | f_lval < 0 | f_uval == Inf | f_uval == -Inf | f_uval < 0) {
    stop("Bad density: density is +/-Inf or negative")
  }
  # check boundness of the bounds on the density
  if (f_lval == 0) {
    stop("Lower bound does not bound the density correctly")
  } else if (f_uval == 0) {
    stop("Upper bound does not bound the density correctly")
  }
  # check lower smaller than upper
  if (lower >= upper) {stop("Lower bound greater than or equal to upper bound")}
  return(TRUE)
}

# test: the three wrong categories
check_support_boundaries(f, 0.2, 0.1)

# check integral of density < Inf
check_density_convergence <- function(f, lower, upper) {
  integral <- try(integrate(f, lower, upper), silent = TRUE)
  if (is(integral, "try-error")) {
    stop("Bad density: density integrates to +/- Inf")
  }
  return(TRUE)
}
check_density_convergence(f, -Inf, 1)


# tests
library(testthat)




