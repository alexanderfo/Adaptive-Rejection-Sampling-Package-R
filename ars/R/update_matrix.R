# update the vertices with new h, h', and secant value
update_vertices <- function(vertices, new_vertex, h){
  new_row <- c(new_vertex, NA, NA, NA)
  vertices <- rbind(vertices, new_row)
  rownames(vertices) <- NULL
  
  # order the vertices matrix by the x-value
  vertices <- vertices[order(vertices[,1]),]
  idx <<- which(vertices[,1] == new_vertex) # global assignment
  
  vertices[idx,2] <- h(new_vertex)
  vertices[idx,3] <- evaluate_deriv(h,new_vertex)
  
  # update secant values
  vertices[idx-1,4] <- calc_secant(vertices, idx-1, idx)
  vertices[idx,4] <- calc_secant(vertices, idx, idx+1)
  
  return(vertices)
}

# update the func_list with u, l, the exponential of u & l, as well as the end points of the u's & l's
update_func_list <- function(vertices, func_list, h, idx){
  # get the lengths of new list
  len_new <- length(vertices[,1])
  len_old <- len_new - 1
  
  # update the u functions
  u_new <- create_u(vertices, idx, h)
  exp_u_new <- exp_fun(u_new)
  func_list$u[idx:len_new] <- c(u_new, func_list$u[idx:len_old])
  func_list$exp_u[idx:len_new] <- c(exp_u_new, func_list$exp_u[idx:len_old])
  
  # update the endpoints/z's of the u functions
  # assuming the new vertex is always between the two starting points
  # that is, idx cannot be 1, or len_old+1
  z_new_left <- calc_intersection(h, vertices, idx-1, idx)
  z_new_right <- calc_intersection(h, vertices, idx, idx+1)
  z_new <- c(z_new_left, z_new_right)
  if(idx == len_old) func_list$z_lo[idx:len_new] <- z_new
  else func_list$z_lo[idx:len_new] <- c(z_new, func_list$z_lo[(idx+1):len_old])
  func_list$z_hi[(idx-1):len_new] <- c(z_new, func_list$z_hi[idx:len_old])
  
  # update the l functions
  l_new_left <- create_l(vertices, idx-1)
  l_new_right <- create_l(vertices, idx)
  func_list$l[(idx-1):len_new] <- c(l_new_left, l_new_right, func_list$l[idx:len_old])
  
  # update the end-points/x's of the l functions
  func_list$x_lo[idx:len_new] <- c(vertices[idx,1], func_list$x_lo[idx:len_old])
  func_list$x_hi[(idx-1):len_new] <- c(vertices[idx,1], func_list$x_hi[(idx-1):len_old])
  
  return(func_list)
}

