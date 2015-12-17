# test
# f: standard normal

# setup
f <- dnorm
h <- function(x) return(log(f(x)))
x <- seq(-5,5,by=0.1)
ymin <- min(h(x))
ymax <- max(h(x))
plot(x, h(x), xlim=c(-5,5), ylim=c(ymin,ymax),type='l',lty=1)

# initialization
vertices <- init_vertices(h, -5, 5)
points(vertices[,1],vertices[,2], pch=20)
u <- update_u(vertices)
l <- update_l(vertices)
lines(x,u(x), lty = 5)
lines(x,l(x), lty = 3)

# for automatic update
update <- function(dat){
  for(pt in dat){
    new_vertex <- pt
    vert_up <- update_vertices(vertices,new_vertex, h)
    vertices <- vert_up$new_vert
    idx <- vert_up$new_idx
    u <- update_u(vertices)
    l <- update_l(vertices)
    
    plot(x, h(x), xlim=c(-5,5), ylim=c(ymin,ymax+3),type='l',lty=1)
    lines(x, u(x), lty=5)
    lines(x, l(x), lty=3)
    points(vertices[,1], h(vertices[,1]), pch = 17)
    len <- length(vertices[,1])
    intersection <- sapply(1:(len-1), function(i) calc_intersection(h, vertices, i, i+1))
    points(intersection, u(intersection), pch = 20)
  }
}

# a few try's
update(0)
update(1)
update(c(-3,2,-0.43,4.2))
update(c(seq(-4.5,4.5,by=1.5)))
update(c(-4.8,-3,-2.5,-2,1,-0.43,0.32,0,2.35))

# for manual update
new_vertex <- 2
vertices <- update_vertices(vertices,new_vertex, h)
u <- update_u(vertices)
l <- update_l(vertices)

plot(x, h(x), xlim=c(-5,5), ylim=c(ymin,ymax+3),type='l',lty=1)
lines(x, u(x), lty=5)
lines(x, l(x), lty=3)
points(vertices[,1], h(vertices[,1]), pch = 17)

z <- get_intersection(h, vertices, -5,5)
  
abline(v=z)
