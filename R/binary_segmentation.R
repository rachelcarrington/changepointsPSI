#' Estimate changepoints in mean, variance, or slope using binary segmentation
#'
#' @param x Vector of data.
#' @param threshold Threshold for determining changepoint candidates; defaults to \code{sqrt(2 * log(length(x)) * mad(x))}.
#' @param maxiter The maximum number of changepoints to find: defaults to the maximum (\code{length(x) - 2} if \code{model = "slope"}, otherwise \code{length(x) - 1}).
#' @param model Changepoint model: should be one of \code{"mean"}, \code{"var"}, or \code{"slope"}. The default is \code{"mean"}.
#' @param loss Loss function to use: for \code{model = "var"} only. One of \code{"cusum"}, \code{"icss"}, or \code{"lrs"} (default). See below for details.
#'
#' @description
#' Binary segmentation algorithm for detecting changes in mean, variance, or slope.
#' The algorithm terminates when either 
#' \enumerate{
#' \item the maximum number of changepoints has been found; or
#' \item there are no remaining points where the CUSUM statistic is above the threshold.
#' }
#'
#' @return A list:
#' \itemize{
#' \item \code{results} Dataframe containing results.
#' \item \code{changepoints} Vector of changepoints detected.
#' \item \code{threshold} Value of \code{threshold} (either specified by user, or calculated in the function).
#' \item \code{maxiter} Value of \code{maxiter} (either specified by user, or default used by function).
#' }
#' 
#' @details
#' There are three options for the changepoint model. In each case the data is assumed to be of the form
#' \eqn{x_t = \mu_t + \sigma_t \epsilon_t}, where \eqn{\epsilon_t \sim N(0, 1)}. The options are:
#' \itemize{
#' \item \code{model = "mean"} (change in mean): \eqn{\mu_t} is constant except at changepoints; \eqn{\sigma_t} is constant
#' \item \code{model = "slope"} (change in slope): \eqn{\mu_t = a_t + t b_t}, with \eqn{a_t} and \eqn{b_t} constant except at changepoints; \eqn{\sigma_t} is constant
#' \item \code{model = "var"} (change in variance): \eqn{\mu_t = 0} and \eqn{\sigma_t} is constant except at changepoints.
#' }
#' For the change in variance model (\code{model = "var"}), there are three possible loss functions for estimating changeoints.
#' If \code{loss = "lrs"}, the likelihood ratio statistic will be used for each interval
#' \deqn{(e - s + 1) \log \sum_{t=s}^e X_t^2 - (\tau - s + 1) \log \sum_{t=s}^{\tau} X_t^2 - (e - \tau) \log \sum_{t=\tau + 1}^e X_t^2}.
#' If \code{loss = "cusum"}, the CUSUM statistic of \eqn{X^2} will be used:
#' \deqn{(\frac{(\tau - s + 1)(e - \tau)}{e - s + 1})^2 (\frac{1}{\tau - s + 1} \sum_{t=s}^{\tau} X_t^2 + \frac{1}{e - \tau} \sum_{t=\tau + 1}^e X_t^2)}
#' If \code{loss = "icss"}, then
#' \deqn{\sum_{t=s}^{\tau} X_t^2 / \sum_{t=s}^e X_t^2 - \frac{\tau}{T}.}
#'
#' @export
#'
#' @examples
#' set.seed(10)
#' x <- c(rep(1,20), rep(-1,30), rep(1,50)) + rnorm(100)
#' binary_segmentation(x)
#'
#' x <- c(seq(1, 10, length.out=50), rep(10, 50)) + rnorm(100)
#' binary_segmentation(x, model="slope")
#'
#' x <- c(rnorm(50, sd=1), rnorm(50, sd=1.5), rnorm(50, sd=2))
#' binary_segmentation(x, model="var", loss="lrs", threshold=5)
#' binary_segmentation(x, model="var", loss="cusum", threshold=10)
#'
binary_segmentation <- function(x, threshold=NULL, maxiter=NULL, model="mean", loss="lrs"){

  stopifnot(is.numeric(x))
  stopifnot(model %in% c("mean", "var", "slope"))
  if ( model == "var" ){
    stopifnot(loss %in% c("lrs", "cusum", "icss"))
  }

  n <- length(x)

  # Set threshold, if unspecified
  if ( is.null(threshold) ){
    if ( model == "mean" ){
      threshold <- sqrt(2 * log(n) * mad(x))
    } else if ( model == "slope" ){
      threshold <- sqrt(2 * log(n) * mad_slope(x))
    } else if ( model == "var" ){
      if ( loss == "cusum" || loss == "cusum0" ){
        threshold <- sqrt(2 * log(n) * max(x^2))
      } else {
        stop("Threshold should be specified if model == var and loss == lrs.")
      }
    }
  } else if ( threshold < 0 ){
    threshold <- 0
  }

  # Set maximum number of iterations, if not specified.
  if ( is.null(maxiter) ){
    if ( model == "mean" || model == "var" ){
      maxiter <- n - 1
    } else if ( model == "slope" ){
      maxiter <- n - 2
    }
  } else if ( maxiter < 1 || maxiter > n - 1 ){
    if ( model == "mean" || model == "var" ){
      warning("Invalid maxiter value, setting to n - 1.")
      maxiter <- n - 1
    } else if ( model == "slope" ){
      warning("Invalid maxiter value, setting to n - 2.")
      maxiter <- n - 2
    }
  }

  results <- matrix(NA, nrow=maxiter, ncol=7)
  colnames(results) <- c("iter", "s", "e", "b", "d", "lrs", "cp")

  # Find first changepoint: calculate the CUSUM/LRS statistics
  iter <- 1
  s <- 1
  e <- n
  s_all <- s
  e_all <- e
  if ( model == "mean" ){
    C0 <- cumsum(x)
    lrs <- cusum(x, cumsums=C0)
  } else if ( model == "slope" ){
#    lambda <- create_lambda_change_in_slope(n, s, e)
#    lambda <- calculate_nu_slope(n, s, e)
#    lrs <- c(0, t(lambda) %*% x)
    lrs <- calculate_lrs_slope(x, s, e)
  } else if ( model == "var" ){
    C0 <- cumsum(x^2)
    if ( loss == "lrs" ){
      lrs <- n * log(C0[n]/n) - (1:(n-1)) * log(C0[-n]/(1:(n-1))) - ((n-1):1) * log((C0[n] - C0[-n])/((n-1):1))
    } else if ( loss == "cusum" ){
      lrs <- cusum(x^2, cumsums=C0)
    } else if ( loss == "icss" ){
      lrs <- cusum(x^2, cumsums=C0, icss=TRUE)
    }
  }
  
  # Find which time point has the maximum value, and check whether the maximum is above the threshold
  b <- which.max(abs(lrs)) # changepoint candidate
  cp <- abs(lrs[b]) > threshold # = 1 if changepoint detected, 0 otherwise
  
  # Find direction of change
  if ( model == "var" & loss == "lrs" ){
    d <- ifelse(C0[b] / C0[n] < b/n, 1, -1)
  } else {
    d <- ifelse(lrs[b] > 0, -1, 1)
  }
  results[1,] <- c(iter, s, e, b, d, lrs[b], cp)

  s_all <- c(s_all, b + 1)
  e_all <- c(e_all, b)

  # Divide at changepoints and repeat
  while ( iter < maxiter ) {
    
    if ( cp == 1 ){
      # If a CP was detected last time
      iter <- iter + 1

      # Re-calculate loss function on (s, e)
      lrs[b] <- NA
      if ( model == "slope" ){
        lrs[b + 1] <- NA
      }
      if ( b > s + 1 ){
        if ( model == "mean" ){
          lrs[s:(b - 1)] <- cusum(x, s=s, e=b, cumsums=C0)
        } else if ( model == "slope" ){
#          lambda <- calculate_nu_slope(n, s, e=b)
#          lrs[(s + 1):(b - 1)] <- t(lambda) %*% x
          lrs[(s + 1):(b - 1)] <- calculate_lrs_slope(x, s, b)[(s + 1):(b - 1)]
        } else if ( model == "var" ){
          if ( loss == "lrs" ){
            if ( s > 1 ){
              lrs[s:(b - 1)] <- (b - s + 1) * log((C0[b] - C0[s - 1])/(b - s + 1)) -
                (1:(b - s)) * log((C0[s:(b - 1)] - C0[s - 1]) / (1:(b - s))) - 
                ((b - s):1) * log((C0[b] - C0[s:(b - 1)])/((b - s):1))
            } else {
              lrs[s:(b - 1)] <- (b - s + 1) * log(C0[b]/b) -
                (1:(b - s)) * log(C0[s:(b - 1)] / (1:(b - s))) - 
                ((b - s):1) * log((C0[b] - C0[s:(b - 1)])/((b - s):1))
            }
          } else if ( loss == "cusum" ){
            lrs[s:(b - 1)] <- cusum(x^2, s=s, e=b, cumsums=C0)
          } else if ( loss == "icss" ){
            lrs[s:(b - 1)] <- cusum(x^2, s=s, e=b, cumsums=C0, icss=TRUE)
          }
        }
      } else if ( b == s + 1 ) {
        lrs[s] <- NA
      }
      if ( e > b + 2 ){
        if ( model == "mean" ){
          lrs[(b + 1):(e - 1)] <- cusum(x, s=b + 1, e=e, cumsums=C0)
        } else if ( model == "slope" ){
#          lambda <- calculate_nu_slope(n, s=b + 1, e=e)
#          lrs[(b + 2):(e - 1)] <- t(lambda) %*% x
          lrs[(b + 2):(e - 1)] <- calculate_lrs_slope(x, b + 1, e)[(b + 2):(e - 1)]
        } else if ( model == "var" ){
          if ( loss == "lrs" ){
            lrs[(b + 1):(e - 1)] <- (e - b) * log((C0[e] - C0[b])/(e - b)) -
              (1:(e - b - 1)) * log((C0[(b + 1):(e - 1)] - C0[b]) / (1:(e - b - 1))) - 
              ((e - b - 1):1) * log((C0[e] - C0[(b + 1):(e - 1)])/((e - b - 1):1))
          } else if ( loss == "cusum" ){
            lrs[(b + 1):(e - 1)] <- cusum(x^2, s=b + 1, e=e, cumsums=C0)
          } else if ( loss == "icss" ){
            lrs[(b + 1):(e - 1)] <- cusum(x^2, s=b + 1, e=e, cumsums=C0, icss=TRUE)
          }
        }
      } else if ( e == b + 2 & e <= n - 1 ){
        lrs[b + 1] <- lrs[b + 2] <- NA
      } else if ( e == b + 1 & e <= n - 1 ){
        lrs[b + 1] <- NA
      } else if ( e == b + 2 & e == n ){
        lrs[b + 1] <- NA
      }

      # Find next changepoint candidate
      b <- which.max(abs(lrs)) 
      cp <- abs(lrs[b]) > threshold
      s <- max(s_all[s_all <= b])
      e <- min(e_all[e_all > b])
      s_all <- c(s_all, b + 1)
      e_all <- c(e_all, b)
      if ( model == "var" & loss == "lrs" ){
        if ( s > 1 ){
          d <- ifelse((C0[b] - C0[s - 1]) / (C0[e] - C0[s - 1]) < b/(e - s + 1), 1, -1)
        } else {
          d <- ifelse(C0[b] / C0[e] < b/e, 1, -1)
        }
      } else {
        d <- ifelse(lrs[b] > 0, -1, 1)
      }
      results[iter,] <- c(iter, s, e, b, d, lrs[b], cp)
        
    } else {
      # If no CP was detected, end while loop
      iter <- maxiter
      results <- results[!is.na(rowSums(results)),,drop=FALSE] # remove empty rows from results matrix
    }
  }

  # Get changepoints
  changepoints <- results[,"b"][results[,"cp"] == 1]

  # Convert results to data frame
  results <- data.frame(results)

  return(list(results=results, changepoints=changepoints, threshold=threshold, maxiter=maxiter))

}