setwd("~/git/stat243-project/")
#setwd("/Users/meikao/Desktop/UC.Berkeley/Academics/STAT243/stat243-project")
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
  
  # check input validity
  check_support_boundaries(density, lb, ub)
  
  h <- function(x, ...) log(density(x, ...))
  mode <- find_mode(density, lb, ub)
  condition <- is_logconcave(h, lb, ub, mode[1], ...)
  #print(condition)
  if(condition == 1) return(runif(n, lb, ub))
  else if(condition == 2) print("Truncated distribution: the leftmost point is the mode.")
  else if(condition == 4) print("Truncated distribution: the right point is the mode.")
  else if(condition == FALSE) stop("Bad density: not log-concave")
  
  # a counter of samples
  numSamples <- 0
  
  # preallocate space for the samples we want to draw
  samples <- rep(0,n) 
  
  # avoid numeric issue and reset the lb and ub if ifinity
  if (condition == 2 && is.infinite(h(ub))) ub = optim(mode[1]+1, function(x) {density(x) - 1e-18}, method = "BFGS")$par
  else {
    if (is.infinite(h(ub))) {
      ub = optim(mode[1]+1, function(x) {density(x) - 1e-18}, method = "BFGS")$par
    } 
    if (is.infinite(h(lb))) {
      lb = optim(mode[1]-1, function(x) {density(x) - 1e-18}, method = "BFGS")$par
    }
  }
  
  # initialize the T_k set in paper
  # vertices <- c(v1, v2), if we decide to use a class
  vertices <- init_vertices(h, lb, ub, condition, mode[1])
  
#   if (is.infinite(lb)) lb = -1000
#   if (is.infinite(ub)) ub = 1000
  u <- update_u(vertices, lb, ub)
  l <- update_l(vertices, h, lb, ub)
  
  while(numSamples < n){
    z <- get_intersection(vertices, lb, ub)
    
    # x, w, bin are vectors
    x <- draw_sample(vertices, z, num_of_samples = n - numSamples)
    # x <- draw_sample(exp_fun(u), z, num_of_samples = n - numSamples)
    w <- runif(n - numSamples)
    
    squeeze <- exp(l(x) - u(x))
    squeeze[is.nan(exp(l(x) - u(x)))] = 0
    
    accept_at_sqz <- (w <= squeeze)
    if(all(accept_at_sqz == TRUE)){
      samples[(numSamples+1):n] <- x
      numSamples <- n
    } 
    else{
      stop_pt <- which(accept_at_sqz == FALSE)[1]
      if(stop_pt > 1){
        diff <- (stop_pt-1) - 1 # the diff is auxiliary so the following 2 slices are of equal length
        samples[(numSamples+1) : (numSamples+1+diff)] <- x[1:(stop_pt - 1)]
        numSamples <- numSamples + (stop_pt - 1)
      }
    }
    
#     if(!check_logconcave_iter(u,l,h,x[stop_pt])){
#       stop(paste("ERROR: The input density is not log-concave at", x[stop_pt]))
#     }
    
    if(numSamples < n){
      accept_at_rej <- (w[stop_pt] <= exp(h(x[stop_pt]) - u(x[stop_pt])))
      if (accept_at_rej) {
        samples[numSamples+1] <- x[stop_pt]
        numSamples <- numSamples + 1
      }
      
      # when the last point finishes the sampling, stop
      if(numSamples < n){
        x_stop_pt <- x[stop_pt]
        update_result <- update_vertices(vertices, x_stop_pt, h)
        vertices <- update_result$vertices
        if (update_result$shrink == TRUE) {
          ifelse(x_stop_pt <= mode[1], lb <- x_stop_pt, ub <- x_stop_pt)
        }
        u <- update_u(vertices, lb, ub)
        l <- update_l(vertices, h, lb, ub)
      }
    }
    if (!all(u(samples[numSamples]) >= l(samples[numSamples])))
      stop("Bad density: not log-concave")
  }
  return(samples)
}
