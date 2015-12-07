# test1
# f: standard normal

# setup
f <- dnorm
h <- function(x) return(log(f(x)))
x <- seq(-5,5,by=0.1)
ymin <- min(h(x))
ymax <- max(h(x))
plot(x, h(x), xlim=c(-7,7), ylim=c(ymin-2,ymax+5),type='l',lty=1)

# initialization
vertices <- init_vertices(h, -5, 5)
points(vertices[,1],vertices[,2], pch=20)
func_list <- init_piecewise(vertices, h, -6, 6)

len <- length(vertices[,1])
for(i in 1:len){
  u <- func_list$u[[i]]
  lines(x, u(x), lty=5)
}
for(i in 1:(len-1)){
  l <- func_list$l[[i]]
  lines(x, l(x), lty=3)
}

# manually update
new_vertex <- 1.96
update_matrix(vertices, func_list, new_vertex, h)
len <- length(vertices[,1])
for(i in 1:len){
  u <- func_list$u[[i]]
  lines(x, u(x), lty=5)
  points(func_list$z_lo[i], func_list$u[[i]](func_list$z_lo[i]), pch = 20)
  points(func_list$z_hi[i], func_list$u[[i]](func_list$z_hi[i]), pch = 19)
}
for(i in 1:(len-1)){
  l <- func_list$l[[i]]
  lines(x, l(x), lty=3)
  points(func_list$x_lo[i], func_list$l[[i]](func_list$x_lo[i]), pch = 20)
  points(func_list$x_hi[i], func_list$l[[i]](func_list$x_hi[i]), pch = 19)
}








