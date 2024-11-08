#' Calculate S for autocorrelated data.
#'
#' @description For a given dataset and changepoint algorithm (one of binary segmentation, wild binary segmentation, or narrowest over threshold),
#' calculate the possible outcomes that can be obtained by applying the algorithm to
#' the perturbed dataset \eqn{X'(\phi)}, and the intervals of \eqn{\phi} values that lead to each possible outcome. See below for more details.
#'
#' @param x Vector of data.
#' @param Sigma Covariance matrix of \code{x}.
#' @param nu Vector that defines the test statistic, which is \eqn{\nu^T x}.
#' @param results Output of changepoint algorithm (either \code{binary_segmentation}, \code{wild_binary_segmentation}, or 
#' \code{narrowest_over_threshold}).
#' @param phi_obs Observed value of the test statistic \eqn{\nu^T x}; optional.
#' @param phi_var Variance of the test statistic; optional.
#' @param first_cp_only Logical. If \code{TRUE}, condition on the fact that the changepoint of interest is in the model; 
#' if \code{FALSE}, condition on all changepoints. Defaults to \code{FALSE}.
#' @param cp_of_interest Index of the changepoint of interest. Defaults to 1.
#'
#' @details
#' \eqn{X'(\phi)} is defined as
#' \deqn{X'(\phi) = X_{obs} + \frac{1}{||\nu||^2} \nu (\phi - \phi_{obs}).}
#'
#' @return A dataframe containing intervals with the changepoints obtained when \eqn{\phi} is in each interval.
#'
#' @export
#'
#' @examples
#' set.seed(100)
#' rho <- 0.7
#' n <- 100
#' x <- arp(n, rho, mu=c(rep(1, n/2), rep(-1, n/2)))
#' results <- find_changepoints(x, "bs")
#' b <- results$changepoints
#' print(b)
#' h <- 10
#' nu <- c(rep(0, b[1] - h), rep(1/h, h), rep(-1/h, h), rep(0, length(x) - b[1] - h))
#' Sigma <- diag((1 - rho^seq(2, 2*n, by=2))/(1 - rho^2))
#' for ( i in 1:(n - 1) ){
#'   Sigma[(i + 1):n, i] <- Sigma[i, (i + 1):n] <- rho^seq(1:(n - i)) * Sigma[i, i]
#' }
#' calculate_S_autocor(x, Sigma, nu, results, first_cp_only=TRUE)
#'
calculate_S_autocor <- function(x, Sigma, nu, results, phi_obs=NULL, phi_var=NULL, first_cp_only=FALSE, cp_of_interest=1){
  
  # Calculate S given x, b, d
  
  n <- length(x)
  
  if ( is.null(phi_obs) ){
    phi_obs <- c(t(nu) %*% x)
  }
  
  b <- results$changepoints
  d <- results$results$d
  maxiter <- results$params$maxiter
  model <- results$params$model
  method <- results$params$method
  cp_params <- within(results$params, rm(model, method))
  
  if ( is.null(phi_var) ){
    phi_var <- c(t(nu) %*% Sigma %*% nu)
  }
  Sigma_nu <- Sigma %*% nu
  
  eps0 <- 0.01
  
  # Calculate limits beyond which we treat intervals as infinity
  upper_limit <- (-1) * qnorm(10^(-20), sd=sqrt(phi_var))
  lower_limit <- qnorm(10^(-20), sd=sqrt(phi_var))
  
  if ( first_cp_only ){
    n.cp <- cp_of_interest
    interval <- calculate_interval(results, nu=nu, Sigma=Sigma, phi_obs=phi_obs, maxiter=n.cp, autocor=TRUE)
    S <- matrix(c(interval, as.matrix(c(b[1:n.cp], d[1:n.cp]))), nrow=1)
    colnames(S) <- c("lower_lim", "upper_lim", paste0("b", 1:n.cp), paste0("d", 1:n.cp))
    
    
    ncp_max <- (ncol(S) - 2)/2
    eps <- eps0
    while ( max(S[,"upper_lim"]) < upper_limit ){
      phi <- max(S[,"upper_lim"]) + eps
      x_phi <- x + (phi - phi_obs) / phi_var * Sigma_nu
      results_phi <- find_changepoints(x_phi, method=method, params=cp_params, model=model)
      b2 <- results_phi$changepoints
      d2 <- results_phi$results$d
      
      if ( b[cp_of_interest] %in% b2 ){ # then we only have to go as far as b[1] in calculating CPs
        n.cp <- (1:length(b2))[b2 == b[cp_of_interest]]
      } else {
        n.cp <- maxiter
      }
      
      interval <- calculate_interval(results_phi, x=x, nu=nu, Sigma=Sigma, phi_obs=phi_obs, phi_var=phi_var, maxiter=n.cp, autocor=TRUE)
      
      # Check this interval is the next one
      ncps_found <- min(n.cp, sum(!is.na(b2)))
      if ( abs(interval[1] - max(S[,"upper_lim"])) < 10^(-10) ){
        if ( ncps_found == 0 ){
          S <- rbind(S, c(interval, rep(NA, ncol(S) - 2)))
        } else if ( ncps_found == ncp_max ){
          S <- rbind(S, c(interval, b2[1:ncps_found], d2[1:ncps_found]))
        } else if ( ncps_found <= ncp_max ){
          S <- rbind(S, c(interval, b2[1:ncps_found], rep(NA, ncp_max - ncps_found), d2[1:ncps_found], rep(NA, ncp_max - ncps_found)))
        } else {
          S <- cbind(S[,1:(2 + ncp_max),drop=FALSE], matrix(NA, nrow=nrow(S), ncol=ncps_found - ncp_max), S[,-(1:(2 + ncp_max)),drop=FALSE],
                      matrix(NA, nrow=nrow(S), ncol=ncps_found - ncp_max))
          S <- rbind(S, c(interval, b2[1:ncps_found], d2[1:ncps_found]))
          ncp_max <- ncps_found
        }
        eps <- eps0
      } else {
        eps <- eps/2
      }
      
    }
    
    eps <- eps0
    while ( min(S[,"lower_lim"]) > lower_limit ){
      phi <- min(S[,"lower_lim"]) - eps
      x_phi <- x + (phi - phi_obs) / phi_var * Sigma_nu
      results_phi <- find_changepoints(x_phi, method=method, params=cp_params, model=model)
      b2 <- results_phi$changepoints
      d2 <- results_phi$results$d
      
      if ( b[cp_of_interest] %in% b2 ){ # then we only have to go as far as b[1] in calculating CPs
        n.cp <- (1:length(b2))[b2 == b[cp_of_interest]]
      } else {
        n.cp <- maxiter
      }
      
      interval <- calculate_interval(results_phi, x=x, nu=nu, Sigma=Sigma, phi_obs=phi_obs, phi_var=phi_var, maxiter=n.cp, autocor=TRUE)
      
      # Check this interval is the next one
      ncps_found <- min(n.cp, sum(!is.na(b2)))
      if ( abs(interval[2] - min(S[,"lower_lim"])) < 10^(-10) ){
        if ( ncps_found == 0 ){
          S <- rbind(S, c(interval, rep(NA, ncol(S) - 2)))
        } else if ( ncps_found == ncp_max ){
          S <- rbind(S, c(interval, b2[1:ncps_found], d2[1:ncps_found]))
        } else if ( ncps_found <= ncp_max ){
          S <- rbind(S, c(interval, b2[1:ncps_found], rep(NA, ncp_max - ncps_found), d2[1:ncps_found], rep(NA, ncp_max - ncps_found)))
        } else {
          S <- cbind(S[,1:(2 + ncp_max),drop=FALSE], matrix(NA, nrow=nrow(S), ncol=ncps_found - ncp_max), S[,-(1:(2 + ncp_max)),drop=FALSE],
                      matrix(NA, nrow=nrow(S), ncol=ncps_found - ncp_max))
          S <- rbind(S, c(interval, b2[1:ncps_found], d2[1:ncps_found]))
          ncp_max <- ncps_found
        }
        eps <- eps0
      } else {
        eps <- eps/2
      }
      
    }
    
    
  } else {
    
    # Calculate interval s.t. the given values of b & d are obtained
    interval <- calculate_interval(results, nu=nu, Sigma=Sigma, phi_obs=phi_obs, phi_var=phi_var, maxiter=maxiter, autocor=TRUE)
    if ( length(b) >= 1 ){
      S <- matrix(c(interval, as.matrix(c(b, d))), nrow=1)
      colnames(S) <- c("lower_lim", "upper_lim", paste0("b", 1:length(b)), paste0("d", 1:length(d)))
    } else {
      S <- matrix(c(interval, as.matrix(c(NA, NA))), nrow=1)
      colnames(S) <- c("lower_lim", "upper_lim", "b1", "d1")
    }
    
    ncp_max <- max(length(b), 1)
    
    # Find other intervals
    eps <- eps0
    while ( max(S[,"upper_lim"]) < upper_limit ){
      phi <- max(S[,"upper_lim"]) + eps
      x_phi <- x + (phi - phi_obs) / phi_var * Sigma_nu
      results_phi <- find_changepoints(x_phi, method=method, params=cp_params, model=model)
      b2 <- results_phi$changepoints
      d2 <- results_phi$results$d
      
      interval <- calculate_interval(results_phi, x=x, nu=nu, Sigma=Sigma, phi_obs=phi_obs, phi_var=phi_var, maxiter=maxiter, autocor=TRUE)
      
      # Check this interval is the next one
      if ( abs(interval[1] - max(S[,"upper_lim"])) < 10^(-10) ){
        if ( is.na(b2[1]) ){
          S <- rbind(S, c(interval, rep(NA, ncol(S) - 2)))
        } else if ( length(b2) == ncp_max ){
          S <- rbind(S, c(interval, b2, d2))
        } else if ( length(b2) <= ncp_max ){
          S <- rbind(S, c(interval, b2, rep(NA, ncp_max - length(b2)), d2, rep(NA, ncp_max - length(b2))))
        } else {
          S <- cbind(S[,1:(2 + ncp_max),drop=FALSE], matrix(NA, nrow=nrow(S), ncol=length(b2) - ncp_max), S[,-(1:(2 + ncp_max)),drop=FALSE],
                      matrix(NA, nrow=nrow(S), ncol=length(b2) - ncp_max))
          S <- rbind(S, c(interval, b2, d2))
          ncp_max <- length(b2)
        }
        eps <- eps0
      } else {
        eps <- eps/2
      }
      
    }
    
    eps <- eps0
    while ( min(S[,"lower_lim"]) > lower_limit ){
      phi <- min(S[,"lower_lim"]) - eps
      x_phi <- x + (phi - phi_obs) / phi_var * Sigma_nu
      results_phi <- find_changepoints(x_phi, method=method, params=cp_params, model=model)
      b2 <- results_phi$changepoints
      d2 <- results_phi$results$d
      
      interval <- calculate_interval(results_phi, x=x, nu=nu, Sigma=Sigma, phi_obs=phi_obs, phi_var=phi_var, maxiter=maxiter, autocor=TRUE)
      
      # Check this interval is the next one
      if ( abs( interval[2] - min(S[,"lower_lim"]) ) < 10^(-10) ){
        if ( is.na(b2[1]) ){
          S <- rbind(S, c(interval, rep(NA, ncol(S) - 2)))
        } else if ( length(b2) == ncp_max ){
          S <- rbind(S, c(interval, b2, d2))
        } else if ( length(b2) <= ncp_max ){
          S <- rbind(S, c(interval, b2, rep(NA, ncp_max - length(b2)), d2, rep(NA, ncp_max - length(b2))))
        } else {
          S <- cbind(S[,1:(2 + ncp_max),drop=FALSE], matrix(NA, nrow=nrow(S), ncol=length(b2) - ncp_max), S[,-(1:(2 + ncp_max)),drop=FALSE],
                      matrix(NA, nrow=nrow(S), ncol=length(b2) - ncp_max))
          S <- rbind(S, c(interval, b2, d2))
          ncp_max <- length(b2)
        }
        eps <- eps0
      } else {
        eps <- eps/2
      }
      
    }
    
  }
  
  S <- data.frame(S)
  colnames(S) <- c("lower_lim", "upper_lim", paste0("b", 1:ncp_max), paste0("d", 1:ncp_max))
  
  return(S)
}