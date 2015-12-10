update_u <- function(vertices){
  len <- length(vertices[,1])
  intersection <- sapply(1:(len-1), function(i) calc_intersection(h, vertices, i, i+1))
  u <- function(x){
    bin_idx <- min(which(intersection >= x))
    if(is.infinite(bin_idx)){
      pt_x <- vertices[len,1]
      pt_y <- vertices[len,2]
      slope <- vertices[len,3]
    }
    else{
      pt_x <- vertices[bin_idx,1]
      pt_y <- vertices[bin_idx,2]
      slope <- vertices[bin_idx,3]
    }
    return(pt_y + slope * (x - pt_x))
  }
  return(Vectorize(u))
}

update_l <- function(vertices){
  l <- function(x){
    bin_idx <- max(which(vertices[,1] <= x))
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
