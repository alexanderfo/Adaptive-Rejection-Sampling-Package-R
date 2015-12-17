#' @title Calculate the intersection of tangent lines
#' 
#' @description Calculate the intersection points of the tangent lines defined by two successive points in the vertices matrix
#'
#' @param vertices The vertices matrix that describe the current vertices
#' @param idx1/idx2 The indices of the two sucessive points that the intersection is calculated for
#' @param eps The tolerance of numerical offset of denominator to avoid invalid devision
#' 
#' @return The abscissa of the intersection points the tangent lines defined by two successive points in the vertices matrix
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

#' @title Get intersection points
#' 
#' @description Calculate and combines intersection points of all tangent lines in the current setup, and append the lower/upper bound of the sampling range
#'
#' @param vertices The vertices matrix that describe the current vertices
#' @param lb/ub The lower/upper bound of the sampling range
#' 
#' @return A numeric vector that stores the absissae of intersection points and lb/ub, in ascending order

get_intersection <- function(vertices, lb, ub){
  len <- length(vertices[,1])
  intersection <- sapply(1:(len-1), function(i) calc_intersection_vert(vertices, i, i+1))
  intersection <- c(lb, intersection, ub)
  return(intersection)
}