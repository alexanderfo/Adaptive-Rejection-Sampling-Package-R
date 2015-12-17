# Initialization, when only upper and lower bound for x is provided and log(density)=log(f(x))=h(x)

init_vertices <- function(h, lb, ub, condition, mode){
  x_lo <- lb
  x_hi <- ub

  row1 <- c(x_lo, h(x_lo), evaluate_deriv(h,x_lo), NA)
  row2 <- c(mode, h(mode), 0 + 1e-8, NA)
  row3 <- c(x_hi, h(x_hi), evaluate_deriv(h,x_hi), NA)

  if (condition == 2) vertices <- rbind(row2, row3)
  else if (condition == 3) vertices <- rbind(row1, row2)
  else vertices <- rbind(row1, row2, row3)
  colnames(vertices) <- c("x", "h(x)", "h_prime(x)", "secant")
  rownames(vertices) <- NULL # remove row names created by rbind

  vertices[1,4] <- calc_secant(vertices, 1, 2)
  if(condition == 4){
    vertices[2,4] <- calc_secant(vertices, 2, 3)
    vertices[3,4] <- vertices[2,4]
  }
  else{
    vertices[2,4] <- vertices[1,4]
  }

  return(vertices)
}
