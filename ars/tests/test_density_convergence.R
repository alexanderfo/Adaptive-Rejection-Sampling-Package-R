
f <- function(x) {return(x^(2-1) * (1-x)^(2-1))}
check_density_convergence(f, -Inf, 1)