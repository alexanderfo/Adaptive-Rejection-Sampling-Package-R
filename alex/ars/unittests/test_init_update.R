setwd("~/git/stat243-project/")
source("ars/R/ars.R")

# test
# f: standard normal

# setup
f <- dnorm
h <- function(x) return(log(f(x)))
x <- seq(-5,5,by=0.1)
ymin <- min(h(x))
ymax <- max(h(x))
plot(x, h(x), xlim=c(-5,5), ylim=c(ymin,ymax),type='l',lty=1)
idx <- 1

# initialization
vertices <- init_vertices(h, -5, 5)
points(vertices[,1],vertices[,2], pch=20)
func_list <- init_piecewise(vertices, h, -6, 10)

len <- length(vertices[,1])
for(i in 1:len){
  u <- func_list$u[[i]]
  lines(x, u(x), lty=5)
}
for(i in 1:(len-1)){
  l <- func_list$l[[i]]
  lines(x, l(x), lty=3)
}

update <- function(dat){
  for(pt in dat){
    # manually update
    new_vertex <- pt
    vert_up <- update_vertices(vertices,new_vertex, h)
    vertices <- vert_up$new_vert
    idx <- vert_up$new_idx
    func_list <- update_func_list(vertices, func_list, h, idx)
    
    len <- length(vertices[,1])
    plot(x, h(x), xlim=c(-5,5), ylim=c(ymin,ymax+3),type='l',lty=1)
    for(i in 1:len){
      u <- func_list$u[[i]]
      dat <- c(func_list$z_lo[i], func_list$z_hi[i])
      lines(dat, u(dat), lty=5)
      points(func_list$z_lo[i], func_list$u[[i]](func_list$z_lo[i]), pch = 20)
      if(i==len) points(func_list$z_hi[i], func_list$u[[i]](func_list$z_hi[i]), pch = 20) 
    }
    
    for(i in 1:(len-1)){
      l <- func_list$l[[i]]
      dat <- c(func_list$x_lo[i], func_list$x_hi[i])
      lines(dat, l(dat), lty=3)
      points(func_list$x_lo[i], func_list$l[[i]](func_list$x_lo[i]), pch = 17)
      if(i==len-1) points(func_list$x_hi[i], func_list$l[[i]](func_list$x_hi[i]), pch = 17)
    }
  }
}

# a few try's
update(0)
update(1)
update(c(seq(-4.5,4.5,by=1.5)))
update(c(-4.8,-3,-2.5,-2,1,-0.43,0.32,0,2.35))

  # manually update
  new_vertex <- -0.5
  vert_up <- update_vertices(vertices,new_vertex, h)
  vertices <- vert_up$new_vert
  idx <- vert_up$new_idx
  func_list <- update_func_list(vertices, func_list, h, idx)
  
  len <- length(vertices[,1])
  plot(x, h(x), xlim=c(-5,5), ylim=c(ymin,ymax),type='l',lty=1)
  for(i in 1:len){
    u <- func_list$u[[i]]
    dat <- c(func_list$z_lo[i], func_list$z_hi[i])
    lines(dat, u(dat), lty=5)
    points(func_list$z_lo[i], func_list$u[[i]](func_list$z_lo[i]), pch = 20)
    if(i==len) points(func_list$z_hi[i], func_list$u[[i]](func_list$z_hi[i]), pch = 20) 
  }
  
  for(i in 1:(len-1)){
    l <- func_list$l[[i]]
    dat <- c(func_list$x_lo[i], func_list$x_hi[i])
    lines(dat, l(dat), lty=3)
    points(func_list$x_lo[i], func_list$l[[i]](func_list$x_lo[i]), pch = 17)
    if(i==len-1) points(func_list$x_hi[i], func_list$l[[i]](func_list$x_hi[i]), pch = 17)
  }
  
