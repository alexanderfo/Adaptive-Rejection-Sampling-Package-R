# Function update matrix, adds new element to vertices matrix and the index for it
# Can be used to then pass it in to update_u, l, and z functions.

update_matrix <- function(vertices, new_vertex,h,h_prime){
  
  new_row<-c(new_vertex)
  
  vertices<-rbind(vertices,new_row)
  
  rownames(vertices)<-NULL
  
  vertices <- vertices[order(vertices[,1]),]
  
  idx<-which(vertices[,1]==new_vertex)
  
  vertices[idx,2]=h
  vertices[idx,3]=h_prime
  
  # find new secant value
  vertices[idx-1,4]=calc_secant(vertices,idx-1,idx)
  vertices[idx,4]=calc_secant(vertices,idx,idx+1)
  
  return(list(matrix=vertices,new_index=idx))
}
