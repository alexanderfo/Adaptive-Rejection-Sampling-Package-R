###### unit tests
test_that("check_support_boundaries works", {
  f <- function(x) {return(x^(2-1) * (1-x)^(2-1))}
  expect_that(check_support_boundaries(f, 0.1, 0.4), is_true())
  expect_that(check_support_boundaries(f, 0.4, 0.1), 
              throws_error("Lower bound greater than or equal to upper bound"))
  expect_that(check_support_boundaries(f, Inf, 0.1), 
              throws_error("Bad bounds: density is +/-Inf or negative at upper/lower bounds"))
  expect_that(check_support_boundaries(f, 0.1, -Inf), 
              throws_error("Bad bounds: density is +/-Inf or negative at upper/lower bounds"))
  expect_that(check_support_boundaries(f, 1.2, 0.1), 
              throws_error("Bad bounds: density is +/-Inf or negative at upper/lower bounds"))
  expect_that(check_support_boundaries(f, -0.1, 0.05), 
              throws_error("Bad bounds: density is +/-Inf or negative at upper/lower bounds"))
  f <- function(x) {return(0)}
  expect_that(check_support_boundaries(f, 0.1, 0.4), 
              throws_error("Bad bounds: density is 0 everywhere within bounds"))
  
})

test_that("check_density_convergence works", {
  f <- function(x) {return(x^(2-1) * (1-x)^(2-1))}
  expect_that(check_density_convergence(f, -Inf, 1),
              throws_error("Bad density: density diverges"))
  expect_that(is.numeric(check_density_convergence(f, 0, 1)), is_true())
})

test_that("find_mode works", {
  expect_that(find_mode(dnorm, -Inf, Inf), equals(0 + sqrt(.Machine$double.eps)))
  expect_that(find_mode(dnorm, -300, 300), equals(0 + sqrt(.Machine$double.eps)))
  expect_that(find_mode(dnorm, -2, -0.5), equals(-0.5 + sqrt(.Machine$double.eps)))
})

