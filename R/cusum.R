#' Calculate CUSUM statistics for a vector.
#'
#' @param x Vector for which to calculate CUSUM statistics.
#' @param s Starting index: defaults to 1.
#' @param e Ending index: defaults to \code{length(x)}.
#' @param icss Whether to use the cumulative sum statistic of Inclan and Tiao (see below), rather than the standard CUSUM statistic. Defaults to \code{FALSE}.
#' @param return_full If \code{TRUE}, the function returns a vector of same length as \code{x}, with indices outside of \code{(s, e - 1)} set to
#' \code{NA}. If \code{FALSE} (default), a vector of length \code{e - s + 1} is returned.
#' @param cumsums Vector of cumulative sums of either \code{x} or \code{x[s:e]}. If not supplied, this will be calculated within the function. 
#' Calculating in advance can save computation time if the function needs to be run multiple times, for example within a changepoint algorithm.
#'
#' @description
#' Calculate CUSUM statistics for a vector \code{x}.
#'
#' @details
#' There are two options for which version of the CUSUM statistic to calculate. If \code{icss = "TRUE"}, the function computes the version used
#' in Iterated Cumulative Sum of Squares (ICSS) -- see the paper by Inclan & Tiao (1994) for details. This is
#' \deqn{C(t) = \sum_{j=1}^t X_j / \sum_{j=1}^T X_j - t/T.} (Note that the paper by Inclan \& Tiao uses the square of the data, in which case this should
#' be supplied to the function as \code{x}.)
#' Otherwise, if \code{icss = "FALSE"}, the function computes the version of the CUSUM statistic used in binary segmentation, etc:
#' \deqn{C(t) = \sqrt{t(T - t)/T} * (1/t \sum_{j=1}^t X_t - 1/(T - t) \sum_{j=t + 1}^T X_t}.
#'
#' @return Vector of CUSUM statistics, of length \code{length(x) - 1} (if \code{return_full = TRUE}), or length \code{e - s} if it is not.
#'
#' @export
#'
#' @examples
#' x <- c(1, 3, 2, 5, 3)
#' cusum(x)
#'
#' @references
#' Carla Inclan & George C. Tiao (1994), Use of cumulative sums of squares for retrospective detection of changes of variance,
#' _Journal of the American Statistical Association_, 89(427), 913--923.
#'
cusum <- function(x, s=1, e=length(x), icss=FALSE, return_full=FALSE, cumsums=NULL){

  stopifnot(s >= 1, s < length(x), e > 1, e <= length(x))

  if ( is.null(cumsums) ){
    cumsums <- cumsum(x[s:e])
  } else if ( length(cumsums) == length(x) ){
    if ( s > 1 ){
      cumsums <- cumsums[s:e] - cumsums[s - 1]
    } else {
      cumsums <- cumsums[s:e]
    }
  } else { # (incorrectly specified cumsums)
    cumsums <- cumsum(x[s:e])
  }

  if ( icss ){
    n <- e - s + 1
    cusum_vec <- cumsums[-n] / cumsums[n] - (1:(n - 1)) / n
  } else {
    n <- e - s + 1
    z <- 1:(n - 1)
    a <- sqrt(z * rev(z) / n)
    cusum_vec <- a * (cumsums[-n] / z - (cumsums[n] - cumsums[-n]) / rev(z))
  }
  
  if ( return_full ){
    cusum_vec <- c(rep(NA, s - 1), cusum_vec, rep(NA, length(x) - e))
  }

  return(cusum_vec)

}
