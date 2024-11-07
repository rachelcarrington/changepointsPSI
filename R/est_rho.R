# Auxiliary function
fun1 <- function(x, k){
  return(1 - x^k)
}

fun2 <- function(x, k){
  return((1 - x^k) * (3 - x^k) / (1 - x^2))
}

#' Estimate autocorrelation parameter and noise variance for AR(1) data.
#'
#' @param x Vector of data
#' @param K Number of lagged variances to calculate. See below.
#' @param rho_grid_size Width of grid size to take proposed \eqn{\rho} values; should be between 0 and 0.5. See below.
#' @param model Changepoint model: either \code{"mean"} (change in mean) or \code{"slope"} (change in slope). Defaults to \code{"mean"}.
#' 
#' @description
#' Estimate the parameters of a model with piecewise constant or piecewise linear mean and AR(1) noise.
#'
#' @details
#' The model for the data is
#' \deqn{X_t = \mu_t + \rho(X_{t-1} - \mu_{t-1}) + \epsilon_t,}
#' where \eqn{\mu_t} is constant except at a small number of changepoints, and \eqn{\epsilon_t \sim N(0, \sigma^2)}.
#' This function estimates \eqn{\rho} and \eqn{\sigma^2} using the method described in Section 4 of Romano et al. (2022). This is based on estimating
#' the variances of the lag-\eqn{k} differenced data (for \eqn{k = 1, \ldots, K}) and minimizing the sum of squared differences
#' between these estimates and the theoretical variances, which are functions of \eqn{\rho} and \eqn{\sigma^2}.
#' 
#' For a given \code{\rho} we can calculate analytically the value of \code{\sigma^2} which minimizes this sum of squares.
#' To estimate both \eqn{\rho} and \eqn{\sigma^2}, the function generates a grid of \eqn{\rho} values between 0 and 1,
#' and calculates the optimal \eqn{\sigma^2} for each. The combination of \eqn{(\rho, \sigma^2)} that gives the overall
#' minimum value to the sum of squares objective is returned as the best estimates for these parameters.
#'
#' The parameter \code{K} is the number of differenced variances to calculate: a larger \code{K} will give a more accurate estimate when
#' the mean of the data is constant, but may lead to biased estimates if there are large or frequent changes within the data, so
#' it is not necessarily optimal to choose a large \code{K}.
#' The parameter \code{rho_grid_size} controls the size of the grid of \eqn{\rho} values: here smaller values should yield more accurate
#' results.
#'
#' @return A list containing the estimated noise variance (\code{sigma2}) and autocorrelation parameter (\code{rho}).
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- arp(500, rho=0.7, sd=2)
#' est_rho(x, 20, 0.01)
#'
#' x <- arp(500, rho=0.7, sd=2, mu=c(rep(1, 200), rep(-1, 200), rep(1, 100)))
#' est_rho(x, 20, 0.01)
#'
#' @references
#' Gaetano Romano, Guillem Rigaill, Vincent Runge and Paul Fearnhead (2022), 
#' Detecting Abrupt Changes in the Presence of Local Fluctuations and Autocorrelated Noise.
#' _Journal of the American Statistical Association_, 117(540), 2147--2162.
#'
est_rho <- function(x, K, rho_grid_size, model="mean"){
  
  if ( model == "mean"){
    
    # Calculate lagged variance estimates
    mads <- sapply(1:K, mad_ar, x=x)
    
    # Generate grid of rho values
    rho_grid <- seq(rho_grid_size, 1 - rho_grid_size, by=rho_grid_size)
    
    # Create vector for sigma values
    sigs <- S_val <- rep(0, length(rho_grid))
    
    # Calculate optimal variance for each rho value
    for ( j in 1:length(rho_grid) ){
      rho_grid_k <- fun1(rho_grid[j], 1:K)
      frac_vals <- rho_grid_k / rho_grid_k[2]
      sigs[j] <- sum(frac_vals * mads) / (2 * sum(frac_vals^2))
      S_val[j] <- sum((2 * sigs[j] * frac_vals - mads)^2)
    }
    
    # Find overall optimum sigma and rho
    sig_est <- sigs[which.min(S_val)]
    rho_hat <- rho_grid[which.min(S_val)]
    
  } else if ( model == "slope" ){
    
    # Calculate lagged variance estimates
    mads <- 6 * sapply(1:K, mad_slope, x=x)
    
    # Generate grid of rho values
    rho_grid <- seq(rho_grid_size, 1 - rho_grid_size, by=rho_grid_size)
    
    # Create vector for sigma values
    sigs <- S_val <- rep(0, length(rho_grid))
    
    # Calculate optimal variance for each rho value
    for ( j in 1:length(rho_grid) ){
      rho_grid_k <- fun2(rho_grid[j], 1:K)
      frac_vals <- rho_grid_k
      sigs[j] <- sum(frac_vals * mads) / (2 * sum(frac_vals^2))
      S_val[j] <- sum((2 * sigs[j] * frac_vals - mads)^2)
    }
    
    # Find overall optimum sigma and rho
    sig_est <- sigs[which.min(S_val)]
    rho_hat <- rho_grid[which.min(S_val)]
    
  }
    
  return(list(sigma2=sig_est, rho=rho_hat))
}