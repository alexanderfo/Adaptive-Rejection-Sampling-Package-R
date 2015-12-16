#setwd("~/git/stat243-project/")
setwd("/Users/meikao/Desktop/UC.Berkeley/Academics/STAT243/stat243-project")
source("ars/R/initialize.R")
source("ars/R/new_draw_sample.R")
source("ars/R/aux_func.R")
source("ars/R/evaluate_deriv.R")
source("ars/R/check_support_boundaries.R")
source("ars/R/check_density_convergence.R")
source("ars/R/new_check_logconcave.R")
source("ars/R/new_update_method.R")
source("ars/R/intersecion.R")
source("ars/R/find_mode.R")

arSampler <- function(density, n, lb = -Inf, ub = Inf, ...){
  print(paste("Sampling ", n, " points now!"))
  
  # check the lower and upper bounds are supportive to the given
  check_support_boundaries(density, lb, ub)
  
  # check the unnormalized density converges
  normcst <- check_density_convergence(density, lb, ub)
  if(!normcst)
    stop("Bad density: density does not converge")
  
  # take the log function
  h <- function(x, ...) log(density(x, ...))
  
  # mode[1] is the abscissa, mode[2] is not used in the major routine
  mode <- find_mode(density, lb, ub)
  
  # consider change the function name, as this does not check concavity any more
  condition <- is_logconcave(h, lb, ub, mode[1], ...)
  if(condition == 1){
    print("Uniform distribution: runif is used to generate sample")
    return(runif(n, lb, ub))
  } 
  else if(condition == 2) print("Truncated distribution: the leftmost point is the mode.")
  else if(condition == 3) print("Truncated distribution: the right point is the mode.")
  
  # a counter of samples
  numSamples <- 0
  # preallocate space for the samples we want to draw
  samples <- rep(0,n) 
  
  # avoid numeric issue and reset the lb and ub if infinity
  # need to include the situation when condition == 3 (when the mode is the ub)
  if (condition == 2 && is.infinite(h(ub))) ub <- optim(mode[1]+1, function(x) {density(x) - 1e-18}, method = "BFGS")$par
  else {
    if (is.infinite(h(ub))) {
      ub <- optim(mode[1]+1, function(x) {density(x) - 1e-18}, method = "BFGS")$par
    } 
    if (is.infinite(h(lb))) {
      lb <- optim(mode[1]-1, function(x) {density(x) - 1e-18}, method = "BFGS")$par
    }
  }
  
  if (!is_logconcave_core(h,lb,ub,TRUE))
    stop("Bad density: not log-concave")
  
  # initialize the vertices set, upper hull, and squeezing functions
  vertices <- init_vertices(h, lb, ub, condition, mode[1])
  u <- update_u(vertices, lb, ub)
  l <- update_l(vertices, h, lb, ub)
  if (!check_local_concave(u, l, vertices[,1]))
    stop("Bad density: not log-concave")
  
  while(numSamples < n){
    z <- get_intersection(vertices, lb, ub)
    
    # x, w, bin are vectors
    x <- draw_sample(vertices, z, num_of_samples = n - numSamples)
    w <- runif(n - numSamples)
    
    squeeze <- exp(l(x) - u(x))
    squeeze[is.nan(squeeze)] <- 0 # avoid numerical issues
    accept_at_sqz <- (w <= squeeze)
    
    if(all(accept_at_sqz == TRUE)){
      samples[(numSamples+1):n] <- x
      numSamples <- n
    }
    else{
      stop_pt <- which(accept_at_sqz == FALSE)[1]
      if(stop_pt > 1){
        diff <- (stop_pt-1) - 1 # the diff is auxiliary so the following 2 slices are of equal length
        samples[(numSamples + 1) : (numSamples + 1 + diff)] <- x[1:(stop_pt - 1)]
        numSamples <- numSamples + (stop_pt - 1)
      }
    }
    
    if(numSamples < n){
      accept_at_rej <- (w[stop_pt] <= exp(h(x[stop_pt]) - u(x[stop_pt])))
      if (accept_at_rej) {
        samples[numSamples+1] <- x[stop_pt]
        numSamples <- numSamples + 1
      }
      
      # if the last point does not finishes the sampling, update the vertices, upper hull, and squeezing function
      # else, stop
      if(numSamples < n){
        x_stop_pt <- x[stop_pt]
        update_result <- update_vertices(vertices, x_stop_pt, h)
        vertices <- update_result$vertices
        if (update_result$shrink == TRUE) {
          if(x_stop_pt <= mode[1]) lb <- x_stop_pt
          else ub <- x_stop_pt
          warning(paste("The user input range is shrinked to: [", lb, ",", ub, "]."))
        }
        u <- update_u(vertices, lb, ub)
        l <- update_l(vertices, h, lb, ub)
        
        if (!check_local_concave(u, l, vertices[,1]))
          stop("Bad density: not log-concave")
      }
    }
#     if (any(u(samples[numSamples]) < l(samples[numSamples])))
#       stop("Bad density: not log-concave")
  }
  return(samples)
}
