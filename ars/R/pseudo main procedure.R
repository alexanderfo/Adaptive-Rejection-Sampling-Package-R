# class definition of a vertex maybe
setClass("vertex",
		representation(
			x = "numeric"
			y = "numeric"
		)
)

arSampler <- function(density, n, lb, ub, ...){
	# check input validity
	# define a modular function validate() to achieve this
	if inputs are invalid{
		print error message
		return()
	} else{
		print "pass" # just to be more interactive, as recommended by Chris
	}

	# a counter of samples
	numSamples <- 0

	# preallocate space for the samples we want to draw
	samples <- rep(0,n) 

	# initialize the T_k set in paper
	# vertices <- c(v1, v2), if we decide to use a class
	vertices <- c((x1, y1), (x2, y2))

	# define modular function blocks to achieve these
	define u_0, s_0, l_0 based on vertices

	while(numSamples <= n){

		x <- drawSample(s, ...)
		w <- runif(1)

		if x satisfy squeeze test{
			accept x
			increment numSamples
		}

		else{
			if x satisfy rejection test{
				accept x
				increment numSamples
			} else{
				reject x
			}

			# define a modular function if we do direct insertion into the right place
			include x into the vertices vector in a sorted manner 
			
			# define modular function blocks to achieve these
			update u(x), s(x), l(x)
		}
	}

	return(samples)
}