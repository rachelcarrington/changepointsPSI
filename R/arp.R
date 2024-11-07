#' Generate data with a piecewise constant mean and AR(p) noise.
#'
#' @param n The number of data points.
#' @param rho A p-dimensional vector of correlation parameters (p >= 1).
#' @param mu The mean of the data -- this can be either a single value or a vector of length \code{n}. Defaults to 0.
#' @param sd The standard deviation of the noise (defaults to 1).
#'
#' @return A vector of length \code{n}.
#' 
#' @details
#' Data is simulated from the model
#' \deqn{X_t = \mu_t + \rho (X_{t-1} - \mu_{t-1}) + \epsilon_t},
#' where \eqn{\epsilon_t \sim N(0, \sigma^2)}.
#' 
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' rho <- 0.6
#' mu <- c(rep(0, 50), rep(1, 50))
#' x <- arp(n, rho, mu)
#' plot(x, type="l")
#'
#' n <- 200
#' rho <- c(0.7, -0.2)
#' mu <- c(rep(0, 50), rep(1, 100), rep(2, 50))
#' x <- arp(n, rho, mu, 0.5)
#' plot(x, type="l")
#' 
arp <- function(n, rho, mu=0, sd=1){

  ord <- length(rho)
  epsilon <- rnorm(n, sd=sd)
  z <- rep(NA, n)
  z[1] <- epsilon[1]
  for ( i in 2:(ord + 1) ){
    z[i] <- sum(rho[1:(i - 1)] * z[(i - 1):1]) + epsilon[i]
  }
  for ( i in (ord + 2):n ){
    z[i] <- sum(rho * z[(i - 1):(i - ord)]) + epsilon[i]
  }
  x <- mu + z
  
  return(x)
}