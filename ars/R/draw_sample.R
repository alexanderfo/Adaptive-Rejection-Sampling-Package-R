draw_sample <- function(envelope, z, num_of_samples = 1){
  # Draw sampling point candidates
  # Args: 
  #   envelope: list of exponentiated linear envelope functions
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