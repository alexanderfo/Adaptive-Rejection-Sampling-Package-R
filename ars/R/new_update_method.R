update_vertices <- function(vertices, new_vertex, h){
  new_row <- c(new_vertex, NA, NA, NA)
  vertices <- rbind(vertices, new_row)
  rownames(vertices) <- NULL
  
  # order the vertices matrix by the x-value
  vertices <- vertices[order(vertices[,1]),]
  idx <- which(vertices[,1] == new_vertex)
  
  vertices[idx,2] <- h(new_vertex)
  vertices[idx,3] <- evaluate_deriv(h,new_vertex)
  
  # update secant values
  vertices[idx-1,4] <- calc_secant(vertices, idx-1, idx)
  vertices[idx,4] <- calc_secant(vertices, idx, idx+1)
  
  # make sure the dummy secant slope of the end point is the same as the second last pt
  len <- length(vertices[,1])
  vertices[len,4] <- vertices[len-1,4]
  
  return(vertices)
}

update_u <- function(vertices, lb, ub){
  intersection <- get_intersection(vertices, lb, ub)
  u <- function(x){
    if(x <= intersection[1] || x >= tail(intersection, 1)){
      return(0)
    }
    else{
      bin_idx <- min(which(intersection > x)) - 1
      pt_x <- vertices[bin_idx,1]
      pt_y <- vertices[bin_idx,2]
      slope <- vertices[bin_idx,3]
    }
    
    if(is.infinite(slope) || is.nan(slope)) return(0)
    return(pt_y + slope * (x - pt_x))
  }
  return(Vectorize(u))
}

update_l <- function(vertices){
  l <- function(x){
    bin_idx <- ifelse(all(vertices[,1] > x), Inf, max(which(vertices[,1] <= x)))
    if(is.infinite(bin_idx) || bin_idx == length(vertices[,1])) return(-Inf)
    else{
      pt_x <- vertices[bin_idx,1]
      pt_y <- vertices[bin_idx,2]
      slope <- vertices[bin_idx,4]
      return(pt_y + slope * (x - pt_x))
    }
  }
  return(Vectorize(l))
}
