# Initialization, when only upper and lower bound for x is provided and log(density)=log(f(x))=h(x)

source("ars/R/evaluate_deriv.R")
source("ars/R/aux_func.R")

init_vertices <- function(h, lb, ub){
#   # defensive programming, part check log-concavity
#   if(is.infinite(x_lo)) {
#     x_lo <- -10^16
#     i_lo <- 16 #exponent for -10^i_lo
#     while(evaluate_deriv(h,x_lo) <= 0) {
#       x_lo <- x_lo - 10^i_lo
#       i_lo <- i_lo + 1
#       if(is.infinite(x_lo)) stop("NOT LOG-CONCAVE")
#     }
#   }
#   
#   if(is.infinite(x_hi)) {
#     x_hi <- 10^16
#     i_hi <- 16
#     while(evaluate_deriv(h,x_lo) >= 0) {
#       x_hi <- x_hi + 10^i_hi
#       i_hi <- i_hi + 1
#       if(is.infinite(x_hi)) stop("NOT LOG-CONCAVE")
#     }
#   }
  
  # assign values to the starting points
  mode <- find_mode(h, lb, ub)[1]
  x_lo <- ifelse(is.infinite(lb), mode - runif(1), lb)
  x_hi <- ifelse(is.infinite(ub), mode - runif(1), ub)
  
  # Build matrix vertices that defines: x values, h(x) values, h_prime(x) value, and the secant slope between x1 to x2 (stored at index 1))
  row1 <- c(x_lo, h(x_lo), evaluate_deriv(h,x_lo), NA)
  row2 <- c(mode, h(mode), evaluate_deriv(h,mode), NA)
  row3 <- c(x_hi, h(x_hi), evaluate_deriv(h,x_hi), NA)
  
  vertices <- rbind(row1, row2, row3)
  colnames(vertices) <- c("x", "h(x)", "h_prime(x)", "secant")
  rownames(vertices) <- NULL # remove row names created by rbind
  
  vertices[1,4] <- calc_secant(vertices, 1, 2)
  vertices[2,4] <- calc_secant(vertices, 2, 3)
  vertices[3,4] <- vertices[2,4]
  
  return(vertices)
}

# input_lo & input_hi are the user-input range of the density
# these are used as z_0 and z_k for initalization
init_piecewise <- function(vertices, h, input_lo, input_hi){
  u1 <- create_u(vertices, 1, h)
  u2 <- create_u(vertices, 2, h)
  l1 <- create_l(vertices, 1)
  
  u_list <- c(u1,u2)
  exp_u_list <- list(exp_fun(u1), exp_fun(u2))
  l_list <- c(l1,NA)
  
  # Initialize z, only creates first value of z, can be made into a function for all other z's
  z1 <- calc_intersection(h, vertices, 1, 2)
  z_lo_vec <- c(input_lo, z1)
  z_hi_vec <- c(z1, input_hi)
  
  x_lo_vec <- c(vertices[1,1], NA)
  x_hi_vec <- c(vertices[2,1], NA)
  
  func_list <- list(u=c(u1,u2), exp_u=exp_u_list, z_lo=z_lo_vec, z_hi=z_hi_vec,
                     l=l_list, x_lo=x_lo_vec, x_hi=x_hi_vec)
  return(func_list)
}
