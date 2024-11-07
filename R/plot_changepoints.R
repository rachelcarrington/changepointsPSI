#' Plot data with estimated changepoint locations and estimated trends.
#'
#' @param results Output of \code{find_changepoints} applied to a vector of data.
#' @param return_mu Whether to return estimates of the trend; defaults to \code{FALSE}.
#'
#' @return If \code{return_mu = TRUE}, vector of estimated parameters (mean, slope, or variance) for each changepoint segment.
#' 
#' @export
#'
#' @examples
#' x <- c(rep(0, 100), rep(1, 100), rep(-1, 100), rep(1, 100), rep(2, 100)) + rnorm(500)
#' results <- find_changepoints(x, "bs")
#' plot_changepoints(results)
#'
plot_changepoints <- function(results, return_mu=FALSE, use_ggplot=FALSE, tvals=NULL){
  
  x <- results$x
  n <- length(x)
  
  if ( is.null(tvals) ){
    tvals <- 1:n
  }
  
  b <- sort(results$changepoints)
  d <- results$results$d
  method <- results$params$method
  model <- results$params$model
  
  if ( model == "mean" ){
    
    if ( use_ggplot ){
      g <- ggplot(data.frame(t=tvals, x=x)) + geom_point(aes(x=t, y=x)) + 
        theme_classic() + labs(x="t", y="X")
    } else {
      plot(tvals, x, xlab="t", ylab="X", pch=16)
    }
    
    if ( length(b) == 0 ){
      mu_vals <- mean(x)
      mu_hat <- rep(mu_vals, n)
      if ( use_ggplot ){
        g <- g + geom_segment(x=tvals[1], xend=tvals[n], y=mu_vals, yend=mu_vals, colour="red")
      } else {
        lines(c(tvals[1], tvals[n]), c(mu_vals, mu_vals), col="red")
      }
    } else {
      mu_vals <- rep(NA, length(b) + 1)
      mu_hat <- rep(NA, n)
      mu_hat[1:b[1]] <- mu_vals[1] <- mean(x[1:b[1]])
      if ( use_ggplot ){
        g <- g + geom_segment(x=tvals[1], xend=tvals[b[1]], y=mu_vals[1], yend=mu_vals[1], colour="red")
      } else {
        lines(c(tvals[1], tvals[b[1]]), c(mu_vals[1], mu_vals[1]), col="red")
      }
      if ( length(b) >= 2 ){
        for ( i in 2:length(b) ){
          mu_hat[(b[i - 1] + 1):b[i]] <- mu_vals[i] <- mean(x[(b[i - 1] + 1):b[i]])
          if ( use_ggplot ){
            g <- g + geom_segment(x=tvals[b[i - 1]], xend=tvals[b[i]], y=mu_vals[i], yend=mu_vals[i], colour="red")
          } else {
            lines(c(tvals[b[i - 1]], tvals[b[i]]), c(mu_vals[i], mu_vals[i]), col="red")
          }
        }
      }
      mu_hat[(b[length(b)] + 1):n] <- mu_vals[length(b) + 1] <- mean(x[(b[length(b)] + 1):n])
      if ( use_ggplot ){
        g <- g + geom_segment(x=tvals[b[length(b)]], xend=tvals[n], y=mu_vals[length(b) + 1], yend=mu_vals[length(b) + 1], colour="red")
      } else {
        lines(c(tvals[b[length(b)]], tvals[n]), c(mu_vals[length(b) + 1], mu_vals[length(b) + 1]), col="red")
      }
    }
    if ( use_ggplot ){
      g <- g + geom_vline(xintercept=results$changepoints, colour="blue", linetype="dashed")
    } else {
      abline(v=results$changepoints, col="blue", lty=2)
    }
    
  } else if ( model == "slope" ){
    
#    C1 <- cumsum(1:n)
#    C2 <- cumsum((1:n)^2)

    b <- sort(results$changepoints)
    K <- length(b)
    M <- matrix(0, K + 2, K + 2)
    X <- rep(0, K + 2)
    for ( i in 1:n ){
      mt <- c(1, i, i - b)
      mt[mt < 0] <- 0
      M <- M + crossprod(t(mt))
      xt_coef <- c(1, i, i - b)
      xt_coef[xt_coef < 0] <- 0
      xt <- xt_coef * x[i]
      X <- X + xt
    }

    theta_hat <- solve(M) %*% X
    
    b <- c(sort(results$changepoints), n)
    mu_hat <- rep(theta_hat[1], n)
    mu_hat <- mu_hat + theta_hat[2] * (1:n)
    for ( i in 2:length(b) ){
      mu_hat[(b[i-1] + 1):n] <- mu_hat[(b[i-1] + 1):n] + theta_hat[i + 1] * (1:(n - b[i-1]))
    }
    if ( use_ggplot ){
      g <- ggplot(data.frame(t=tvals, x=x, mu=mu_hat)) + geom_point(aes(x=t, y=x)) + 
           geom_line(aes(x=t, y=mu), colour="red") +
           geom_vline(xintercept=results$changepoints, col="blue", linetype="dashed") +
           theme_classic() + labs(x="t", y="X")
    } else {
      plot(tvals, x, xlab="t", ylab="X")
      lines(tvals, mu_hat, col="red")
      abline(v=results$changepoints, col="blue", lty=2)
    }
    
  } else if ( model == "var" ){
    
    plot(1:n, x, xlab="t", ylab="X", pch=16)
    abline(v=results$changepoints, col="blue", lty=2)
    
  }
  
  if ( !use_ggplot ){
    g <- NULL
  }
  
  if ( return_mu ){
    if ( model == "mean" ){
      return(list(mu_vals, g))
    } else if ( model == "slope" ){
      return(list(theta_hat, g))
    } else if ( model == "var" ){
      b_sorted <- c(0, sort(results$changepoints), n)
      var_vals <- rep(NA, length(b_sorted) - 1)
      for ( i in 1:(length(b_sorted) - 1) ){
        var_vals[i] <- var(x[(b_sorted[i] + 1):b_sorted[i + 1]])
      }
      return(list(var_vals, g))
    }
  } else {
    return(list(NULL, g))
  }
  
}
