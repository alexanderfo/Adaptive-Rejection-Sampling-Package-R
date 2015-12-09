arSampler <- function(density, n, lb = -Inf, ub = Inf){
	# check input validity
	# define a modular function validate() to achieve this
	#if inputs are invalid{
	#	print error message
	#	return()
	#} else{
	#	print "pass" # just to be more interactive, as recommended by Chris
	#}

	# a counter of samples
	numSamples <- 0

	# preallocate space for the samples we want to draw
	samples <- rep(0,n) 
	
	# a "global" variable keeping track of the most recent vertex
	idx <- 1

	# initialize the T_k set in paper
	# vertices <- c(v1, v2), if we decide to use a class
	vertices <- init_vertices(density, lb, ub)
	# define modular function blocks to achieve these
	func_list <- init_piecewise(vertices, density, lb, ub)

	while(numSamples <= n){
	  len <- func_list$z_lo
	  intersection <- c(func_list$z_lo, func_list$z_hi[len])
	  
	  # x, w, bin are vectors
		x <- draw_sample(func_list$exp_u, intersection, num_of_samples = n - numSamples)
		w <- runif(n - numSamples)
    bin <- find_bin(x, func_list$z_hi)
    
    # change to lapplyyyyyyy
    for(i in 1:length(bin)){
      l_i <- func_list$l[[bin[i]]]
      u_i <- func_list$u[[bin[i]]]
      squeeze[i] <- exp(l_i(x[i]) - u_i(x[i]))
    }
		
    accept_at_sqz <- w < squeeze
    if(all(accept_at_sqz == TRUE)){
      samples[(numSamples+1):n] <- x
      numSamples <- n
    } 
    else{
      stop_pt <- which(accept_at_sqz == FALSE)[1]
      if(stop_pt > 1){
        diff <- (stop_pt-1) - 1
        samples[(numSamples+1):(numSamples+1+diff)] <- x[1:(stop_pt - 1)]
        numSamples <- numSamples + (stop_pt - 1)
      }
    }
    
    if(numSamples < n){
      if (w[stop_pt] <= exp(h(x[stop_pt]) - func_list$u[[bin[stop_pt]]])) {
        samples[(numSamples+1)] <- x[stop_pt]
        numSamples <- numSamples + 1
      }
      
      # when the last point finishes the sampling, stop
      if(numSamples < n){
        vertices <- update_vertices(vertices, x[stop_pt], h)
        vertices <- update_func_list(vertices, func_list, h, idx)
      }
    }
	}

	return(samples)
}