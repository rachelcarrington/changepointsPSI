kernel_exponential <- function(x1, x2, l=1){
  # Given two vectors x1 and x2, calculate K(x1, x2)
  dists <- abs(outer(x1, x2, "-"))
  return(exp(-0.5 / l^2 * dists))
}