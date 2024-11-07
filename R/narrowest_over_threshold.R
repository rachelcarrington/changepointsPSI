

#' Narrowest over threshold changepoint algorithm.
#'
#' @description
#' Apply the narrowest over threshold changepoint algorithm to a vector of data. Can be used to detect changes in mean, slope, or variance.
#'
#' @param x Vector of data.
#' @param num_rand_samples Number of random intervals to use (defaults to 1000). Ignored if \code{random_samples} is specified.
#' @param random_samples N x 2 matrix containing random intervals (optional; can be generated within the function).
#' @param threshold Changepoint detection threshold, should be a non-negative number. Defaults to \code{sqrt(2 * log(length(x)) * mad(x))}.
#' @param maxiter Maximum number of changepoints to return: defaults to \code{n - 1} (or \code{n - 2} if \code{model = "slope"}.
#' @param model Changepoint model: one of \code{"mean"} (change in mean), \code{"slope"} (change in slope), or \code{"var"} (change in variance).
#' @param loss Loss function (only if \code{model = "var"}): one of \code{"lrs"}, \code{"cusum"}, or \code{"icss"}. (See below.)
#'
#' @return A list.
#' \itemize{
#' \item \code{results} Dataframe containing results
#' \item \code{changepoints} Vector of changepoints detected
#' \item \code{random_samples} N x 2 matrix containing random intervals used in the algorithm.
#' \item \code{threshold} Value of \code{threshold}
#' \item \code{maxiter} Value of \code{maxiter}
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
#' set.seed(1)
#' x <- rnorm(200) + c(rep(1,90), rep(-1,20), rep(1,90))
#' results <- narrowest_over_threshold(x, num_rand_samples=500, threshold=4)
#' print(results$results)
#' print(results$changepoints)
#'
narrowest_over_threshold <- function(x, num_rand_samples=1000, random_samples=NULL, threshold=NULL, maxiter=NULL, model="mean", loss="lrs", return_full_results=FALSE){

  stopifnot(is.numeric(x))
  stopifnot(model %in% c("mean", "var", "slope"))
  if ( model == "var" ){
    stopifnot(loss %in% c("lrs", "cusum", "icss"))
  }

  n <- length(x)

  # Set maximum number of iterations (defaults to n - 1, if unspecified)
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

  # Set threshold, if unspecified
  if ( is.null(threshold) ){
    threshold <- sqrt(2 * log(n) * mad(x))
  } else if ( threshold < 0 ){
    threshold <- 0
  }

  # Draw random samples
  if ( is.null(random_samples) ){
    if ( model == "mean" || model == "var" ){
      random_samples <- generate_random_intervals(n, num_rand_samples, 2)
    } else if ( model == "slope" ){
      random_samples <- generate_random_intervals(n, num_rand_samples, 3)
    }
  } else {
    num_rand_samples <- nrow(random_samples)
  }

  # Sort intervals by width
  interval_widths <- random_samples[,2] - random_samples[,1] + 1
  random_samples <- random_samples[order(interval_widths), ]
  interval_widths <- sort(interval_widths)
  random_samples <- random_samples[interval_widths >= ifelse(model == "slope", 3, 2),, drop=FALSE] # remove any intervals which are too narrow

  results_full <- matrix(NA, nrow=nrow(random_samples), ncol=7)
  colnames(results_full) <- c("index", "s", "e", "b", "d", "lrs", "width")

  # For each interval, calculate CUSUM/LRS statistic
  if ( model == "var" ){
    C0 <- cumsum(x^2)
  } else if ( model == "mean" ){
    C0 <- cumsum(x)
  }
  for ( interval_index in 1:nrow(random_samples) ){
    s <- random_samples[interval_index, 1]
    e <- random_samples[interval_index, 2]
    if ( model == "mean" ){
      lrs <- cusum(x, s, e, cumsums=C0)
      b <- which.max(abs(lrs))
      d <- ifelse(lrs[b] < 0, 1, -1)
      results_full[interval_index,] <- c(interval_index, s, e, b + s - 1, d, lrs[b], (e - s + 1))
    } else if ( model == "slope" ){
#      lambda <- calculate_nu_slope(n, s, e, return_full=FALSE)
#      lrs <- t(lambda) %*% x[s:e]
      lrs <- calculate_lrs_slope(x, s, e, return_full=FALSE)
      b <- which.max(abs(lrs))
      d <- ifelse(lrs[b] < 0, 1, -1)
      results_full[interval_index,] <- c(interval_index, s, e, b + s, d, lrs[b], (e - s + 1))
    } else if ( model == "var" ){
      lrs <- rep(0, n - 1)
      if ( loss == "lrs" ){
        if ( s > 1 ){
          lrs[s:(e - 1)] <- (e - s + 1) * log((C0[e] - C0[s - 1])/(e - s + 1)) - (1:(e - s)) * log((C0[s:(e - 1)] - C0[s - 1]) / (1:(e - s))) -
            ((e - s):1) * log((C0[e] - C0[s:(e - 1)])/((e - s):1))
        } else {
          lrs[s:(e - 1)] <- (e - s + 1) * log(C0[e]/e) - (1:(e - s)) * log(C0[s:(e - 1)] / (1:(e - s))) - 
            ((e - s):1) * log((C0[e] - C0[s:(e - 1)])/((e - s):1))
        }
      } else if ( loss == "cusum" ){
        lrs[s:(e - 1)] <- cusum(x^2, s, e, cumsums=C0)
      } else if ( loss == "icss" ){
        lrs[s:(e - 1)] <- cusum(x^2, s, e, cumsums=C0, icss=TRUE)
      }
      b <- which.max(abs(lrs))
      if ( loss == "lrs" ){
        if ( s > 1 ){
          d <- ifelse((C0[b] - C0[s - 1])/(b - s) < (C0[e] - C0[s - 1])/(e - s + 1), 1, -1)
        } else {
          d <- ifelse(C0[b] / b < C0[e] / e, 1, -1)
        }
      } else {
        d <- ifelse(lrs[b] < 0, 1, -1)
      }
      results_full[interval_index,] <- c(interval_index, s, e, b, d, lrs[b], (e - s + 1))
    }
  }
  changepoint_candidates <- results_full[abs(results_full[,"lrs"]) > threshold, , drop=FALSE]
  

  if ( length(changepoint_candidates) > 0 ){
  
    # Within each width, sort by CUSUM statistic
    widths <- unique(changepoint_candidates[,"width"])
    for ( width in widths ){
      z2 <- changepoint_candidates[changepoint_candidates[,"width"] == width,,drop=FALSE]
      if ( nrow(z2) > 1 ){
        changepoint_candidates[changepoint_candidates[,"width"] == width, ] <- z2[order(abs(z2[,"lrs"]), decreasing=TRUE), ]
      }
    }
    
    # Create empty dataframe
    results <- numeric(0)
    cps_found <- 0

    # Get results
    while( length(changepoint_candidates) > 0 & cps_found < maxiter ){
      cps_found <- cps_found + 1
      results <- rbind(results, changepoint_candidates[1,])
      changepoint_candidates <- changepoint_candidates[!(changepoint_candidates[,"s"] <= changepoint_candidates[1, "b"] & changepoint_candidates[,"e"] > changepoint_candidates[1,"b"]),,drop=FALSE]
    }
    
    # Convert to dataframe
    results <- data.frame(results)
    colnames(results) <- c("index", "s", "e", "b", "d", "lrs", "width")
  } else {
    # If no changepoints were detected, create empty dataframe
    results <- data.frame(index=NA, s=NA, e=NA, b=NA, d=NA, lrs=NA, width=NA)
  }
  
  if ( !return_full_results ){
    results_full <- NA
  }

  return( list(results=results, changepoints=results$b, random_samples=random_samples, results_full=results_full, threshold=threshold, maxiter=maxiter) )
}
