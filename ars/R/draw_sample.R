draw_sample <- function(vertices, z, num_of_samples = 1) {
  # Draw sampling point candidates
  # Args: 
  #   vertices: matrix of x, h(x), h'(x), secant(x) (n points)
  #   z: vector of all intersection points including lower and upper bounds (n+1)
  # Output:
  #   candidates: vector sampling candidates
  numerator <- exp(vertices[, 2]) * (exp((z[-1] - vertices[, 1]) * vertices[, 3]) 
                            - exp((z[-length(z)] - vertices[, 1]) * vertices[, 3])) / vertices[, 3]
  denominator <- sum(numerator)
  cdf.zs <- c(0, numerator / denominator)
  w <- runif(num_of_samples)
  indices <- sapply(w, function(w_i) {sum(cumsum(cdf.zs) < w_i)})
  candidates <- (log(numerator[indices] * runif(num_of_samples) * vertices[indices, 3] + exp(vertices[indices, 2] + (z[indices] - vertices[indices, 1]) * vertices[indices, 3])) 
                 - vertices[indices, 2]) / vertices[indices, 3] + vertices[indices, 1]
  return(candidates)
}