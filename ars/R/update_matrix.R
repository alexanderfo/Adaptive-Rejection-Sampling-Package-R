# Function that updates the vertices matrix and return the row index of the updated element
# The updated matrix can be used to pass in to the function update_u, update_l and update_z functions

update_matrices <- function(vertices, func_list, new_vertex, h){
  
  # update the vertices
  new_row<-c(new_vertex, NA, NA, NA)
  
  vertices<-rbind(vertices,new_row)
  
  rownames(vertices)<-NULL
  
  vertices <- vertices[order(vertices[,1]),]
  
  idx<-which(vertices[,1]==new_vertex)
  
  vertices[idx,2]=h(new_vertex)
  vertices[idx,3]=evaluate_deriv(h,new_vertex)
  
  # find new secant value
  vertices[idx-1,4]=calc_secant(vertices,idx-1,idx)
  vertices[idx,4]=calc_secant(vertices,idx,idx+1)
  
  # update the df
  len <- length(vertices[,1])
  u_new <- create_u(vertices, idx, h)
  exp_u_new <- exp_fun(u_new)
  func_list$u[idx:len] <- c(u_new, func_list$u[idx:(len-1)])
  func_list$exp_u[idx:len] <- c(u_new, func_list$exp_u[idx:(len-1)])
  
  z_new_left <- calc_intersection(h, vertices, idx-1, idx)
  z_new_right <- calc_intersection(h, vertices, idx, idx+1)
  
  z_lo_replace <- c(z_new_left, z_new_right, func_list$z_lo[(idx+1):(len-1)])
  truncate_idx <- which(is.na(z_lo_replace))
  z_lo_replace <- z_lo_replace[1:truncate_idx-1]
  func_list$z_lo[idx:len] <- z_lo_replace
  func_list$z_hi[(idx-1):len] <- c(z_new_left, z_new_right, func_list$z_hi[idx:(len-1)])
  
  l_new_left <- create_l(vertices, idx-1)
  l_new_right <- create_l(vertices, idx)
  
  l_replace <- c(l_new_left, l_new_right)
  if(idx+1 < len) l_replace <- c(l_replace, func_list$l[(idx):(len-1)]) # if the new point is at the last
  func_list$l[(idx-1):len] <- l_replace
  
  x_lo_replace <- c(new_vertex, func_list$x_lo[idx:(len-1)])
  x_hi_replace <- c(new_vertex, func_list$x_hi[(idx-1):(len-1)])
  func_list$x_lo[idx:len] <- x_lo_replace
  func_list$x_hi[(idx-1):len] <- x_hi_replace
  
  return(list(new_vertices=vertices, new_df=func_list, new_index=idx))
}
