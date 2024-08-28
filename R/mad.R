#' Estimate variance of a vector by median absolute deviation, assuming piecewise constant mean
#'
#' @description Estimate variance of a vector using median absolute deviation, under the assumption that the data 
#' is generated with a piecewise constant mean plus Gaussian noise.
#'
#' @param x Vector of data.
#' @param lag Lag at which to calculate differences; defaults to 1.
#'
#' @return Estimated variance of \code{x}.
#' @export
#'
#' @examples
#' x <- rnorm(200, sd=2) + c(rep(1, 100), rep(-1, 100))
#' mad(x)
#'
mad <- function(x, lag=1){

  difs <- x[-(1:lag)] - x[1:(length(x) - lag)]
  mad_est <- median(abs(difs - median(difs))) / qnorm(0.75)
  mad_est <- mad_est^2 / 2

  return(mad_est)
}


#' Estimate variance of a vector by median absolute deviation, assuming piecewise linear mean
#'
#' @description Estimate variance of a vector using median absolute deviation, under the assumption that the data 
#' is generated with a piecewise linear mean plus Gaussian noise.
#'
#' @param x Vector of data.
#' @param lag Lag at which to calculate differences; defaults to 1.
#'
#' @return Estimated variance of \code{x}.
#' @export
#'
#' @examples
#' x <- rnorm(200, sd=2) + seq(0, 3, length.out=200)
#' mad_slope(x)
#'
mad_slope <- function(x, lag=1){

  difs <- x[-(1:(2 * lag))] - 2 * x[-c(1:lag, (length(x) - lag + 1):length(x))] + x[-((length(x) - 2 * lag + 1):length(x))]
  mad_est <- median(abs(difs)) / qnorm(0.75)
  mad_est <- mad_est^2 / 6

  return(mad_est)
}
