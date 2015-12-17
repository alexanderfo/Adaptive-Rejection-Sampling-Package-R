#' @title Draw sampling point candidates
#' 
#' @description Draw samples from the density described by the upper hull function
#'
#' @param vertices The vertices matrix that describe the current vertices
#' @param z The vector of abscissa of the intersection points of tangent lines
#' @param num_of_samples The number of samples to draw
#' 
#' @return The candidate samples that are drawn from the upper hull functions, which will soon be tested and accepted/rejected
draw_sample <- function(vertices, z, num_of_samples = 1) {
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