update_vertices <- function(vertices, new_vertex, h){
  # At new_vertex, density != 0 AND derivative exists, o/w shrink bounds and
  # update shrink vertices accordingly
  h_new_vertex <- h(new_vertex)
  hp_new_vertex <- evaluate_deriv(h,new_vertex)
  
  if (is.infinite(h_new_vertex) || is.infinite(hp_new_vertex) || is.nan(hp_new_vertex)) {
    return(list(vertices = vertices, shrink = TRUE))
  }
  
  # otherwise, update the vertices
  new_row <- c(new_vertex, NA, NA, NA)
  vertices <- rbind(vertices, new_row)
  rownames(vertices) <- NULL
  
  # order the vertices matrix by the x-value
  vertices <- vertices[order(vertices[,1]),]
  idx <- which(vertices[,1] == new_vertex)
  
  vertices[idx,2] <- h_new_vertex
  vertices[idx,3] <- hp_new_vertex
  
  # update secant values
  len <- length(vertices[,1])
  if(idx > 1) vertices[idx-1,4] <- calc_secant(vertices, idx-1, idx)
  if(idx < len) vertices[idx,4] <- calc_secant(vertices, idx, idx+1)
  
  # make sure the dummy secant slope of the end point is the same as the second last pt
  # this happens only when the last secant value is changed (idx >= len - 1)
  if(idx >= len-1) vertices[len,4] <- vertices[len-1,4]
  
  return(list(vertices = vertices, shrink = FALSE))
}

update_u <- function(vertices, lb, ub){
  intersection <- get_intersection(vertices, lb, ub)
  u <- function(x){
    # when x is out of bounds
    lb_new <- intersection[1]
    ub_new <- intersection[length(intersection)]
    if(x < lb_new || x > ub_new) return(0)
    
    # determine which tangent line should x belong to
    if(abs(x - ub_new) < sqrt(.Machine$double.eps)) bin_idx <- length(intersection) - 1 # avoid nan value from which() fucnction
    else bin_idx <- min(which(intersection > x)) - 1
    
    # find the x, y and slope values of the point that determines the tangent line
    pt_x <- vertices[bin_idx,1]
    pt_y <- vertices[bin_idx,2]
    slope <- vertices[bin_idx,3]
    
    if(is.infinite(slope) || is.nan(slope)) return(0)
    return(pt_y + slope * (x - pt_x))
  }
  return(Vectorize(u))
}

update_l_old <- function(vertices, h, lb, ub){
  l <- function(x){
    # when x is out of bounds of vertices
    #if(x < vertices[1,1] || x > vertices[length(vertices[,1]),1]) return(-Inf)
    num_vertices <- nrow(vertices)
    if (x < vertices[1,1] ) {
      pt_x <- lb
      ifelse(is.infinite(h(lb)), pt_y <- -1e8, pt_y <- h(lb))
      slope <- (vertices[1,2] - pt_y) / (vertices[1,1] - lb)
      return(pt_y + slope * (x - pt_x))
    } else if (x > vertices[num_vertices,1]) {
      pt_x <- ub
      ifelse(is.infinite(h(ub)), pt_y <- -1e8, pt_y <- h(ub))
      slope <- (pt_y - vertices[num_vertices,2]) / (ub - vertices[num_vertices,1])
      return(vertices[num_vertices,2] + slope * (x - vertices[num_vertices,1]))
    } else {
      bin_idx <- max(which(vertices[,1] <= x))
      pt_x <- vertices[bin_idx,1]
      pt_y <- vertices[bin_idx,2]
      slope <- vertices[bin_idx,4]
      return(pt_y + slope * (x - pt_x))
    }
  }
  return(Vectorize(l))
}

update_l <- function(vertices, h, lb, ub){
  l <- function(x){
    # when x is out of bounds of vertices
    if(x < vertices[1,1] || x > vertices[length(vertices[,1]),1]) return(-Inf)
    else {
      bin_idx <- max(which(vertices[,1] <= x))
      pt_x <- vertices[bin_idx,1]
      pt_y <- vertices[bin_idx,2]
      slope <- vertices[bin_idx,4]
      return(pt_y + slope * (x - pt_x))
    }
  }
  return(Vectorize(l))
}
