# calculate the z values, i.e. the intersection points
# strict log-concavity ensures that the denominator is never 0 analytically
# non-strict will be treat in the wrapper
calc_intersection_vert <- function(vertices, idx1, idx2, eps = 1e-8){
  x_lo <- vertices[idx1,1]
  x_hi <- vertices[idx2,1]
  numerator <- vertices[idx2,2] - vertices[idx1,2] - x_hi * vertices[idx2,3] + x_lo * vertices[idx1,3]
  denominator <- vertices[idx1,3] - vertices[idx2,3]
  
  if(is.infinite(vertices[idx1,3]) || is.nan(vertices[idx1,3])) return(x_lo)
  else if(is.infinite(vertices[idx2,3]) || is.nan(vertices[idx2,3])) return(x_hi)
  else if(denominator < eps) return(x_lo)
  else return(numerator/denominator)
}

get_intersection <- function(vertices, lb, ub){
  len <- length(vertices[,1])
  intersection <- sapply(1:(len-1), function(i) calc_intersection_vert(vertices, i, i+1))
#   if(is.infinite(lb)) lb <- -38
#   if(is.infinite(ub)) ub <- 38
  intersection <- c(lb, intersection, ub)
  return(intersection)
}