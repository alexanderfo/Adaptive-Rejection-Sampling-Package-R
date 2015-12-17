#' @title Adaptive Rejection Sampler
#' 
#' @description Draw n samples from a log-concave density using adaptive rejection sampling
#'
#' @param density A log-concave density function (could be unnormlalized)
#' @param n Number of samples
#' @param lb,ub The lower / upper bound of the domain, where the points are sampled from
#' 
#' @return A vector of values that are sampled from the density input within the lb and ub.
ars <- function(density, n, lb = -Inf, ub = Inf, ...){
  # check the lower and upper bounds are supportive to the given
  check_support_boundaries(density, lb, ub)
  
  # update mode, finite lb, finite ub
  init_pts <- find_mode(density, lb, ub)
  mode <- init_pts[1]
  lb <- init_pts[3]
  ub <- init_pts[4]
  
  # check the density is not 0 everywhere
  f_lval <- density(lb)
  f_uval <- density(ub)
  if (is.infinite(log(f_lval)) && is.infinite(log(f_uval)) &&
      is.infinite(log(init_pts[2] - f_lval))) {
    stop("Bad bounds: density is 0 everywhere within bounds")
  }
  
  # check the unnormalized density converges
  normcst <- check_density_convergence(density, lb, ub)
  if(!normcst)
    stop("Bad density: density does not converge")
  
  # take the log function
  h <- function(x, ...) log(density(x, ...))
  
  # determine the shape of the density
  condition <- is_logconcave_shape(h, lb, ub, mode[1], ...)
  if(condition == 1){
    return(runif(n, lb, ub))
  } 
  
  # a counter of samples
  numSamples <- 0
  # preallocate space for the samples we want to draw
  samples <- rep(0,n) 
  
  # avoid numeric issue and reset the lb and ub if infinity
  if (is.infinite(h(ub))) {
    ub <- optim(mode[1]+1, function(x) {density(x) - 1e-18}, method = "BFGS")$par
  } 
  if (is.infinite(h(lb))) {
    lb <- optim(mode[1]-1, function(x) {density(x) - 1e-18}, method = "BFGS")$par
  }
  
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
  }
  return(samples)
}