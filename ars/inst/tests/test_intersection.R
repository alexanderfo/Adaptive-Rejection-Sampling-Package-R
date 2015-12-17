# test

# setup 1
# dnorm
f <- dnorm
h <- function(x) return(log(f(x)))
lb <- -5
ub <- 5
x <- seq(lb-1, ub+1, by=0.1)
ymin <- min(h(x))
ymax <- max(h(x))
plot(x, h(x), xlim=c(lb-1,ub+1), ylim=c(ymin,ymax),type='l',lty=1)
points(c(lb,ub), h(c(lb,ub)), pch=20)

vertices <- init_vertices(h, -5, 5, 3)
z <- get_intersection(vertices, -5, 5)
u <- update_u(vertices, -5, 5)
l <- update_l(vertices)

lines(x, u(x), lty=5)
lines(x, l(x), lty=3)

vertices <- update_vertices(vertices, 1.96, h)
plot(x, h(x), xlim=c(lb-1,ub+1), ylim=c(ymin,ymax+5),type='l',lty=1)
points(vertices[,1], h(vertices[,1]), pch=20)
u <- update_u(vertices, -5, 5)
l <- update_l(vertices)
lines(x, u(x), lty=5)
lines(x, l(x), lty=3)

z <- get_intersection(vertices, lb, ub)

