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
  # check input validity
  check_support_boundaries(density, lb, ub)
  
  h <- function(x, ...) log(density(x, ...))
  mode <- find_mode(density, lb, ub)
  condition <- is_logconcave(h, lb, ub, mode[1], ...)
  print(condition)
  if(condition == 1) return(runif(n, lb, ub))
  else if(condition == 2) print("do something")
  else if(condition == FALSE) stop("Bad density: not log-concave")
  
  # a counter of samples
  numSamples <- 0
  
  # preallocate space for the samples we want to draw
  samples <- rep(0,n) 
  
  # initialize the T_k set in paper
  # vertices <- c(v1, v2), if we decide to use a class
  vertices <- init_vertices(h, lb, ub, condition)
  u <- update_u(vertices, lb, ub)
  l <- update_l(vertices)
  
  while(numSamples < n){
    z <- get_intersection(vertices, lb, ub)
    
    # x, w, bin are vectors
    x <- draw_sample(exp_fun(u), z, num_of_samples = n - numSamples)
    w <- runif(n - numSamples)
    
    squeeze <- rep(0, length(x))
    
    # sapply alternative of the for loop is here
    # please rewrite it if there is smarter ways
    #squeeze <- sapply(1:length(x), 
    #                  function(i){
    #                    l_i <- func_list$l[[l_bin[i]]]
    #                    u_i <- func_list$u[[u_bin[i]]]
    #                    return(exp(l_i(x[i]) - u_i(x[i])))
    #                    })
    
    for(i in 1:length(x)){
      squeeze[i] <- exp(l(x[i]) - u(x[i]))
    }
    
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
        vertices <- update_vertices(vertices, x[stop_pt], h)
        u <- update_u(vertices, lb, ub)
        l <- update_l(vertices)
      }
    }
  }
  
  return(samples)
}
