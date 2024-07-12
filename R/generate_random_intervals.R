#' Generate random intervals.
#'
#' @description Generate random intervals within the range \code{[1, n]}. For use in wild binary segmentation and narrowest over threshold algorithms.
#'
#' @details
#' If the number of intervals supplied is greater than the total number of possible unique intervals contained in \eqn{[1, n]} with width at least
#' \code{min_width}, then the output will be a list of all possible intervals contained in \eqn{[1, n]}.
#'
#' @param n Upper limit.
#' @param N Number of intervals to generate.
#' @param min_width Minimum width of interval allowed (defaults to 2).
#'
#' @return N x 2 matrix of integers
#' @export
#'
#' @examples
#' generate_random_intervals(100, 10, 2)
#'
generate_random_intervals <- function(n, N, min_width=2){

  stopifnot(is.numeric(n) || is.integer(n))
  stopifnot(is.numeric(N) || is.integer(N))
  stopifnot(is.numeric(min_width) || is.integer(min_width))
  stopifnot(min_width >= 0, min_width <= n)

  # Maximum number of random intervals for given n:
  N_max <- (n - min_width + 1) * (n - min_width + 2) / 2

  if ( N < N_max ){

    rand_ints <- matrix(sample(1:n, 2 * N, replace=TRUE), ncol=2)

    # Delete and replace duplicates
    rand_ints <- t(apply(rand_ints, 1, sort))
    rand_ints <- rand_ints[!(duplicated(rand_ints)), ]
    while( nrow(rand_ints) < N ){
      rand_ints_2 <- matrix(sample(1:n, 2*(N - nrow(rand_ints)), replace=TRUE), ncol=2)
      rand_ints_2 <- t(apply(rand_ints_2, 1, sort))
      rand_ints <- rbind(rand_ints, rand_ints_2)
      rand_ints <- rand_ints[!(duplicated(rand_ints)), ]
    }

    # Delete and replace any intervals which are smaller than min_width
    for ( i in 1:N ){
      if ( abs(rand_ints[i,1] - rand_ints[i,2]) < (min_width - 1) ){
        if ( rand_ints[i,1] <= (min_width - 1) ){
          rand_ints[i,2] <- sample((rand_ints[i,1] + min_width - 1):n, 1)
        } else if ( rand_ints[i,1] >= (n - min_width + 1) ){
          rand_ints[i,2] <- sample(1:(rand_ints[i,1] - min_width + 1), 1)
        } else {
          rand_ints[i,2] <- sample(c(1:(rand_ints[i,1] - min_width + 1), (rand_ints[i,1] + min_width - 1):n), 1)
        }
      }

      # Ensure that for each interval the lower limit is in column 1, upper limit in column 2
      if ( rand_ints[i,1] > rand_ints[i,2] ){
        rand_ints[i,] <- rand_ints[i,2:1]
      }
    }

  } else {

    # Return all possible intervals
    ri2 <- matrix(rep(min_width:n, n - min_width + 1), nrow=n - min_width + 1)
    rand_ints <- cbind(rep(1:(n - min_width + 1), (n - min_width + 1):1), as.vector(ri2[lower.tri(ri2, diag=TRUE)]))

  }

  return(rand_ints)

}
