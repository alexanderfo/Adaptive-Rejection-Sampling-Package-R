setwd("~/git/stat243-project/")
source("ars/R/new_update_method.R")
source("ars/R/intersecion.R")

# test

# setup 1
# dnorm
f <- dnorm
h <- function(x) return(log(f(x)))
x <- seq(-5, 5, by=0.1)
ymin <- min(h(x))
ymax <- max(h(x))
plot(x, h(x), xlim=c(-5,5), ylim=c(-10,ymax),type='l',lty=1)

# setup 1
# dnorm
f <- dnorm
h <- function(x) return(log(f(x)))
x <- seq(-5, 5, by=0.1)
ymin <- min(h(x))
ymax <- max(h(x))
plot(x, h(x), xlim=c(-5,5), ylim=c(-10,ymax),type='l',lty=1)
# initialization
vertices <- init_vertices(h, -5, 5)
points(vertices[,1],vertices[,2], pch=20)
u <- update_u(vertices)
l <- update_l(vertices)
lines(x,u(x), lty = 5)
lines(x,l(x), lty = 3)

# setup 2
# exp legit
f <- dexp
h <- function(x) return(log(f(x)))
x <- seq(0.5, 5, by = 0.1)
ymin <- min(h(x))
ymax <- max(h(x))
plot(x, h(x), xlim=c(0,5), ylim=c(-6,ymax),type='l',lty=1)

# initialization
vertices <- init_vertices(h, 0.5, 5)
points(vertices[,1],vertices[,2], pch=20)
u <- update_u(vertices)
l <- update_l(vertices)
lines(x,u(x), lty = 5)
lines(x,l(x), lty = 3)

# for manual update
new_vertex <- 3
vert_up <- update_vertices(vertices,new_vertex, h)
vertices <- vert_up$new_vert
idx <- vert_up$new_idx
u <- update_u(vertices)
l <- update_l(vertices)

plot(x, h(x), xlim=c(0.5,5), ylim=c(ymin,ymax+3),type='l',lty=1)
lines(x, u(x), lty=5)
lines(x, l(x), lty=3)
points(vertices[,1], h(vertices[,1]), pch = 17)

z <- get_intersection(h,vertices, 5)
abline(v=z, lty=3)

# setup 3
# laplace
f <- laplace_pdf
h <- function(x) return(log(f(x)))
x <- seq(-5, 5, by = 0.1)
ymin <- min(h(x))
ymax <- max(h(x))
plot(x, h(x), xlim=c(-5,5), ylim=c(-6,ymax),type='l',lty=1)

# initialization
vertices <- init_vertices(h, -4, 5)
points(vertices[,1],vertices[,2], pch=20)
u <- update_u(vertices)
l <- update_l(vertices)
lines(x,u(x), lty = 5)
lines(x,l(x), lty = 3)

# for manual update
new_vertex <- 3
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

points(intersection, u(intersection), pch = 20)

z <- get_intersection(h,vertices, 5)
abline(v=z, lty=3)
