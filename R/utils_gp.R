# Functions for implementing Gaussian process

kernel_exponential <- function(x1, x2, l=1){
  # Given two vectors x1 and x2, calculate K(x1, x2)
  dists <- abs(outer(x1, x2, "-"))
  return(exp(-0.5 / l^2 * dists))
}

calculate_posterior <- function(x, x0, y0, l=1){
  # Calculate estimate of the posterior distribution of x
  
  # (x0, y0) are observed values
  # x are additional x values, as many as possible
  # l is a hyperparameter
  
  K00 <- kernel_exponential(x0, x0, l)
  K0 <- kernel_exponential(x0, x, l)
  K <- kernel_exponential(x, x, l)
  K00_inv <- solve(K00)
  
  mu_out <- t(K0) %*% K00_inv %*% y0
  cov_out <- K - t(K0) %*% K00_inv %*% K0
  
  return(list(mu=mu_out, cov=cov_out))
}