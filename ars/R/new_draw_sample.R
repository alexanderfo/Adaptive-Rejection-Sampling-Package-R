draw_sample <- function(u, z, num_of_samples = 1){
  # Draw sampling point candidates
  # Args: 
  #   u: function of upper hull (n)
  #   z: vector of all intersection points including lower and upper bounds (n+1)
  # Output:
  #   candidates: vector sampling candidates
  v_integrate <- Vectorize(integrate, vectorize.args = c("lower", "upper"))
  numerator <- v_integrate(u, z[1:(length(z)-1)], z[-1])
  numerator <- unlist(numerator[1, ])
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
    integrate(u, lower = z[w_i], upper = x)$value/denominator - w[i] + ifelse(w_i == 1, 0, cdf_end_pts[w_i-1])
  }
  candidates <- sapply(1:num_of_samples, 
                       function(i) {
                         uniroot(f,interval = c(z[w_idx[i]], z[w_idx[i]+1]), 
                                 i, w_idx[i], tol = 1e-6)$root})
  return(candidates)
}