#' Wild Binary Segmentation
#'
#' @description Implements wild or seeded binary segmentation algorithm for changepoint detection.
#'
#' @param x Vector of data.
#' @param num_rand_samples num_rand_samples Number of random intervals to use (defaults to 1000). Ignored if \code{random_samples} is specified.
#' @param random_samples N x 2 matrix of random intervals to use (optional; can be generated within the function).
#' @param threshold Changepoint detection threshold, should be a non-negative number. Defaults to \code{sqrt(2 * log(length(x)) * mad(x))}.
#' @param maxiter Maximum number of changepoints to find: defaults to \code{n - 1} (or \code{n - 2} if \code{model = "slope"}.
#' @param model Changepoint model: one of \code{"mean"} (change in mean; default), \code{"slope"} (change in slope), or \code{"var"} (change in variance).
#' @param loss Loss function (only if \code{model = "var"}): one of \code{"lrs"}, \code{"cusum"}, or \code{"icss"}. (See below.)
#' @param seeded Whether to implement seeded binary segmentation (defaults to \code{FALSE}).
#' @param decay Decay parameter for seeded binary segmentation; only used if \code{seeded = TRUE}.
#' @param return_full_results If \code{TRUE}, a data frame containing the changepoint candidate for each random interval will be returned with the results.
#' @param prev_full_results Data frame or matrix containing the full output from a previous run of the algorithm.
#' This can be used to save computation time if the algorithm is run repeatedly, with the same set of random intervals, on perturbations of the data; otherwise, it should be left \code{NULL}.
#' @param hlims Vector of length 2, which defines the boundaries of the region of interest.
#' Used if the algorithm is run repeatedly on perturbations of the original data, where data are perturbed within a window \code{(hlims[1], hlims[2])}.
#' This means that the values of the CUSUM statistic corresponding to random intervals entirely outside of this window do not change, and so do not have to be re-calculated.
#'
#' @return A list.
#' \itemize{
#' \item \code{results} Dataframe containing results.
#' \item \code{changepoints} Vector of changepoints detected.
#' \item \code{rand_ints} N x 2 matrix containing random intervals used in the algorithm.
#' \item \code{threshold} Value of \code{threshold} used.
#' \item \code{maxiter} Value of \code{maxiter} used.
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
#' \deqn{(e - s + 1) \log \sum_{t=s}^e X_t^2 - (\tau - s + 1) \log \sum_{t=s}^{\tau} X_t^2 - (e - \tau) \log \sum_{t=\tau + 1}^e X_t^2.}
#' If \code{loss = "cusum"}, the CUSUM statistic of \eqn{X^2} will be used:
#' \deqn{(\frac{(\tau - s + 1)(e - \tau)}{e - s + 1})^2 (\frac{1}{\tau - s + 1} \sum_{t=s}^{\tau} X_t^2 + \frac{1}{e - \tau} \sum_{t=\tau + 1}^e X_t^2).}
#' If \code{loss = "icss"}, then
#' \deqn{\sum_{t=s}^{\tau} X_t^2 / \sum_{t=s}^e X_t^2 - \frac{\tau}{T}.}
#' The parameters \code{return_full_results}, \code{prev_full_results}, and \code{hlims} are used inside the function \code{calculate_pvals_all}, but can otherwise be ignored.
#'
#' @export
#'
#' @examples
#' set.seed(100)
#' x <- rnorm(100) + c(rep(1,40), rep(0,10), rep(1,50))
#' results <- wild_binary_segmentation(x, num_rand_samples=500, threshold=3)
#' print(results$results)
#'
wild_binary_segmentation <- function(x, num_rand_samples=1000, random_samples=NULL, threshold=NULL, maxiter=NULL, model="mean", loss="lrs",
                                        seeded=FALSE, decay=sqrt(2), return_full_results=FALSE, prev_full_results=NULL, hlims=NULL){

  stopifnot(is.numeric(x))
  stopifnot(model %in% c("mean", "var", "slope"))
  if ( model == "var" ){
    stopifnot(loss %in% c("lrs", "cusum", "icss"))
  }

  n <- length(x)

  # Set threshold, if unspecified
  if ( is.null(threshold) ){
    threshold <- sqrt(2 * mad(x) * log(n))
  } else if ( threshold < 0 ){
    threshold <- 0
  }

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

  # Draw random/seeded samples
  if ( is.null(random_samples) ){
    if ( seeded ){
      random_samples <- seeded.intervals(n, decay)
      num_rand_samples <- nrow(random_samples)
    } else {
      if ( is.null(num_rand_samples) ){
        num_rand_samples <- 1000
      }
      if ( model == "mean" || model == "var" ){
        random_samples <- generate_random_intervals(n, num_rand_samples, 2)
      } else if ( model == "slope" ){
        random_samples <- generate_random_intervals(n, num_rand_samples, 3)
      }
    }
  } else {
    num_rand_samples <- nrow(random_samples)
  }

  # Calculate changepoint estimate for each interval
  results_full <- matrix(NA, nrow=num_rand_samples, ncol=7)
  colnames(results_full) <- c("index", "s", "e", "b", "d", "lrs", "cp")
  inds <- 1:num_rand_samples
  if ( !is.null(prev_full_results) & !is.null(hlims) ){
    results_full <- prev_full_results
    
    # For random intervals outside the window of interest, keep previous results
    inds_to_keep <- c(inds[random_samples[inds, 2] < hlims[1]], inds[random_samples[inds, 1] > hlims[2]])
    results_full[inds_to_keep,] <- prev_full_results[inds_to_keep,]
    inds <- inds[!(inds %in% inds_to_keep)]
#    inds <- inds[!(random_samples[inds, 2] < hlims[1])]
#    inds <- inds[!(random_samples[inds, 1] > hlims[2])]
#    results_full <- prev_full_results[!(1:num_rand_samples %in% inds),]
  }
  
  if ( model == "var" ){
    C0 <- cumsum(x^2)
  } else if ( model == "mean" ){
    C0 <- cumsum(x)
  }

  for ( ind in inds ){
    s <- random_samples[ind, 1]
    e <- random_samples[ind, 2]
    if ( model == "mean" ){
      lrs <- cusum(x, s, e, cumsums=C0)
    } else if ( model == "slope" ){
      lambda <- calculate_nu_slope(n, s, e, return_full=FALSE)
      lrs <- t(lambda) %*% x[s:e]
    } else if ( model == "var" ){
      if ( loss == "lrs" ){
        if ( s > 1 ){
          lrs <- (e - s + 1) * log((C0[e] - C0[s - 1])/(e - s + 1)) - (1:(e - s)) * log((C0[s:(e - 1)] - C0[s - 1]) / (1:(e - s))) -
            ((e - s):1) * log((C0[e] - C0[s:(e - 1)])/((e - s):1))
        } else {
          lrs <- (e - s + 1) * log(C0[e]/e) - (1:(e - s)) * log(C0[s:(e - 1)] / (1:(e - s))) - 
            ((e - s):1) * log((C0[e] - C0[s:(e - 1)])/((e - s):1))
        }
      } else if ( loss == "cusum" ){
        lrs <- cusum(x^2, s, e, cumsums=C0)
      } else if ( loss == "icss" ){
        lrs <- cusum(x^2, s, e, cumsums=C0, icss=TRUE)
      }
    }
    b <- which.max(abs(lrs))
    cp <- abs(lrs[b]) > threshold
    if ( model == "var" & loss == "lrs" ){
      if ( s > 1 ){
        d <- ifelse((C0[s + b - 1] - C0[s - 1])/b < (C0[e] - C0[s - 1])/(e - s + 1), 1, -1)
      } else {
        d <- ifelse(C0[b] / b < C0[e] / e, 1, -1)
      }
    } else {
      d <- ifelse(lrs[b] > 0, -1, 1)
    }
    results_full[ind, ] <- c(ind, s, e, b + s - 1, d, lrs[b], cp)
  }

#  if ( return_full_results ){
#    results0_orig <- results_full[order(results_full[,"index"]),]
#  } else {
#    results0_orig <- NA
#  }
  
  # Drop intervals where the CUSUM statistic is below the threshold
  changepoint_candidates <- results_full[results_full[,"cp"] == 1,,drop=FALSE]
  
  if( nrow(changepoint_candidates) == 0 ){

    results <- data.frame(index=NA, s=NA, e=NA, b=NA, d=NA, lrs=NA, cp=NA)
    
  } else {
    # Order by absolute value of test statistics
   changepoint_candidates <- changepoint_candidates[order(abs(changepoint_candidates[,"lrs"]), decreasing=TRUE),,drop=FALSE]
    
    results <- numeric(0)
    iter <- 1
    while( iter <= maxiter ){
      # Take the changepoint from the remaining interval with the highest LRS value
      results <- rbind(results, changepoint_candidates[1,])
      b <- changepoint_candidates[1, "b"]
      cp <- changepoint_candidates[1, "cp"]
      
      # Remove intervals which contain this changepoint
      contains_changepoint <- (b >= changepoint_candidates[,"s"]) & (b < changepoint_candidates[,"e"])
      changepoint_candidates <- changepoint_candidates[!(contains_changepoint),,drop=FALSE]
      
      # Check whether there are any intervals left; if not, end loop
      if ( is.null(changepoint_candidates) ){
        iter <- maxiter + 1
      } else if ( nrow(changepoint_candidates) == 0 ){
        iter <- maxiter + 1
      } else {
        iter <- iter + 1
      }
    }
    
    results <- data.frame(results)
  }
  
  if ( return_full_results ){
    return(list(results=results, changepoints=results$b, random_samples=random_samples, threshold=threshold, maxiter=maxiter, results_full=results_full))
  } else {
    return(list(results=results, changepoints=results$b, random_samples=random_samples, threshold=threshold, maxiter=maxiter, results_full=NA))
  }
}
