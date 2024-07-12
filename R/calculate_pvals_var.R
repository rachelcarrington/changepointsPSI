#' Calculate estimates of p-values for changepoint methods, using Gaussian processes and/or importance sampling.
#'
#' @param results Output of \code{find_changepoints}.
#' @param h Window size: the null hypothesis is that there are no changepoints within a window \code{h} of the changepoint of interest. This can be a positive integer (>=2) 
#' or a vector of length 2 (if the window is asymmetric) or \code{NULL} (if the window is defined by the changepoints either side of the changepoint of interest).
#' @param N_sample Number of samples to use in estimating the p-value; defaults to 100.
#' @param frac Proportion of samples to be used to fit the Gaussian process; the remainder will be used for importance sampling from the Gaussian process posterior.
#' If \code{frac = 1} (which is the default), only a Gaussian process will be used. If \code{frac = 0}, only importance #' sampling will be used; 
#' the sampling will be from either a Beta or Uniform distribution (according to the parameter \code{sampling_distribution}).
#' @param sampling_distribution Distribution from which to sample values of \eqn{\phi} to fit the Gaussian process; or, if \code{frac = 0}, the distribution from which to 
#' sample \eqn{\phi} for importance sampling. Options are \code{"uniform"} (default) or \code{"beta"}. 
#' @param scaling_parameter Scaling parameter for the sampling distribution, if \code{sampling_distribution = "beta"}. The sampling is
#' done from \eqn{Beta(h/2k, h/2k)}. Defaults to 1. Ignored if \code{sampling_distribution = "uniform"}.
#' @param cp_bound If \code{TRUE}, then if there is an estimated changepoint within the window that is fixed under the null hypothesis,
#' this changepoint will be used as the boundary of the window.
#' @param gamma Value in \code{(0, 1]}. If \code{h = NULL}, then \eqn{h} is taken to be \code{gamma} times the distance between the changepoint
#' of interest and the estimated changepoints on either side. Ignored if \code{h} specified; defaults to \code{1}.
#' @param num_pvals Number of p-values to calculate; by default p-values will be calculated for all detected changepoints.
#'
#' @return A list.
#' \itemize{
#' \item \code{changepoints} A vector of detected changepoints.
#' \item \code{p} A vector of estimated p-values.
#' \item \code{P_phi_in_S} A matrix containing estimates of \eqn{Pr(\phi \in S)} for each changepoint and each sampled \eqn{W}.
#' \item \code{p_W} A matrix containing values of the p-value for each changepoint and each sampled \eqn{W}.
#' }
#'
#' @details
#' This function estimates p-values for the change in variance model, when changepoint locations are estimated using methods based on the likelihood ratio: 
#' in this case p-values cannot be calculated analytically. If changepoints are estimated using the CUSUM statistic, then the function \code{calculate_pvals_all}
#' should be used to calculate exact p-values (although this function can still be used to obtain estimates).
#' The function \code{find_changepoints} should first be run on the data to estimate changepoint locations, and the output of that function fed to \code{calculate_pvals_var}.
#'
#' Given a changepoint of interest \eqn{\tau_j}, there are several options for the null hypothesis:
#' \itemize{
#' \item There are no changepoints within a window size \eqn{h} of \eqn{\tau_j}. In this case \code{h} should be supplied. It is possible
#' to use different values of \code{h} on either side of \eqn{\tau_j}, in which case \code{h} should be a vector of length 2. If \eqn{\tau_j - h < 1},
#' then the value of \code{h} to the left of \eqn{\tau_j} will be set to \eqn{\tau_j}; similarly, if \eqn{\tau_h + h > n}, then the value of \code{h}
#' to the right of \eqn{\tau_j} will be set to \eqn{n - \tau_j}.
#' If \code{cp_bound = TRUE}, then if there is another estimated changepoint within \code{h} of \eqn{\tau_j} then the window size used
#' will be the minimum of \code{h} and the distance between \eqn{\tau_j} and the closest estimated changepoint to \eqn{\tau_j}.
#' \item There are no other changepoints between \eqn{\tau_{j-1}} and \eqn{\tau_{j+1}}. In this case \code{h} should be set to \code{NULL}
#' and \code{gamma} should be 1. (This is the default for the function.)
#' \item There are no other changepoints between the midpoint of \eqn{\tau_{j-1}} and \eqn{\tau_j} and the midpoint of \eqn{\tau_j} and 
#' \eqn{\tau_{j+1}}. In this case use \code{h = NULL} and \code{gamma = 0.5}.
#' }
#'
#' @export
#'
#' @examples
#' set.seed(12)
#' x <- c(rnorm(100), rnorm(100, sd=2), rnorm(100))
#' results_wbs <- find_changepoints(x, method="wbs", model="var", 
#'                                  params=list(num_rand_samples=200, threshold=7))
#' pvals <- calculate_pvals_var(results_wbs, h=50, N_sample=200, frac=0.8)
#' plot(x)
#' abline(v=pvals$changepoints, col=ifelse(pvals$p < 0.05, 2, 1))
#' 
#' results_pelt <- find_changepoints(x, method="pelt", model="var", 
#'                                   params=list(penalty="Manual", pen.value=7))
#' pvals <- calculate_pvals_var(results_pelt, h=50, N_sample=200, frac=0.8)
#' plot(x)
#' abline(v=pvals$changepoints, col=ifelse(pvals$p < 0.05, 2, 1))
#' 
calculate_pvals_var <- function(results, h, N_sample=100, frac=1, sampling_distribution="uniform", scaling_parameter=1, NW=1, cp_bound=TRUE, gamma=1, num_pvals=NA){

  # ! need to update with our method?
  
  x0 <- results$x
  n <- length(x0)
  model <- results$params$model
  method <- results$params$method
  params <- within(results$params, rm(model, method))
  b <- results$changepoints
  
  if ( is.na(num_pvals) || num_pvals > length(b) ){
    num_pvals <- length(b)
  }
  
  if ( frac < 0 ){
    frac <- 0
  } else if ( frac > 1 ){
    frac <- 1
  }
  
  N_for_GP <- round(N_sample * frac)
  N_for_IS <- N_sample - N_for_GP
#  p_hat <- rep(NA, num_pvals)
  P_phi_in_S <- p_hat <- matrix(NA, nrow=num_pvals, ncol=NW)
  
  for ( cp_index in 1:num_pvals ){
    
    # Calculate h1 and h2
    hvals <- get_h1_h2(h, b, n, b[cp_index], cp_bound, gamma)
    h1 <- hvals[1]
    h2 <- hvals[2]
    nh <- h1 + h2
    
    # Calculate phi_obs (the observed value of the test statistic) and the critical values of phi
    phi_obs <- sum(x0[(b[cp_index] - h1 + 1):b[cp_index]]^2) / sum(x0[(b[cp_index] - h1 + 1):(b[cp_index] + h2)]^2)
    phi_star <- qbeta(1 - pbeta(phi_obs, h1/2, h2/2), h1/2, h2/2)
    phi_upper <- max(phi_obs, phi_star)
    phi_lower <- min(phi_obs, phi_star)
    
    # Generate psi (W) values
    Psi <- get_psi_vals(x0, h1, h2, NW, b[cp_index], include_original=TRUE, model="var", C0=sum(x[(b[cp_index] - h1 + 1):(b[cp_index] + h2)]^2), 
                        C1=sum(x[(b[cp_index] - h1 + 1):b[cp_index]]^2))$Psi
    
    for ( psi_index in 1:NW ){
      
      x <- x0
      x[(b[cp_index] - h1 + 1):(b[cp_index] + h2)] <- Psi[psi_index,]
      
      # Train Gaussian process
      if ( N_for_GP > 0 ){
        
        # Sample phi values for fitting GP
        if ( sampling_distribution == "uniform" ){
          phi <- runif(N_for_GP) / N_for_GP + seq(0, (N_for_GP - 1)/N_for_GP, length.out=N_for_GP)
        } else if ( sampling_distribution == "beta" ){
          z <- runif(N_for_GP) / N_for_GP + seq(0, (N_for_GP - 1)/N_for_GP, length.out=N_for_GP)
          phi <- qbeta(z, h1/(2 * scaling_parameter), h2/(2 * scaling_parameter))
        }
        
        # Calculate I(phi in S) for each phi
        alpha <- sqrt(phi/phi_obs)
        beta <- sqrt((1 - phi)/(1 - phi_obs))
        inS <- rep(NA, N_for_GP)
        for ( j in 1:N_for_GP ){
          x_phi <- x
          x_phi[(b[cp_index] - h1 + 1):(b[cp_index] + h2)] <- c(alpha[j] * x[(b[cp_index] - h1 + 1):b[cp_index]], beta[j] * x[(b[cp_index] + 1):(b[cp_index] + h2)])
          results_phi <- find_changepoints(x_phi, method=method, params=params, model=model)
          if ( length(results_phi$changepoints) >= 1 ){
            inS[j] <- ifelse(b[cp_index] %in% results_phi$changepoints, 1, 0)
          } else {
            inS[j] <- 0
          }
        }
        
        # If no phi samples are in S, return NA, as the p-value
        if ( sum(inS) == 0 ){
          
          warning("No phi samples were in S, try running with a larger N.")
          p_hat[cp_index, psi_index] <- NA
#          P_both[cp_index, psi_index] <- NA
          P_phi_in_S[cp_index, psi_index] <- NA
          
        # otherwise, keep going 
        } else {
          
          # Fit GP
          X <- sort(c(phi, runif(1000)))
          posterior_results <- calculate_posterior(X, phi, inS, l=100)
          mu_pos <- posterior_results$mu
          mu_pos[mu_pos < 0] <- 0
          mu_pos[mu_pos > 1] <- 1
          
          # If we are not also doing IS, use the GP estimate of the p-value
          if ( N_for_IS == 0 ){
            
            pi_hat <- mu_pos * dbeta(X, h/2, h/2)
#            P_both[cp_index, psi_index] <- sum(pi_hat[X <= phi_lower]) + sum(pi_hat[X >= phi_upper])
            P_phi_in_S[cp_index, psi_index] <- mean(pi_hat)
            p_hat[cp_index, psi_index] <- (sum(pi_hat[X <= phi_lower]) + sum(pi_hat[X >= phi_upper])) / sum(pi_hat)
#            P_phi_in_S[cp_index, psi_index] <- mean(mu_pos)
            
          # Importance sampling (post-GP)
          } else {
            
            # Calculate f/g ratio for our initial phi values
            if ( sampling_distribution == "uniform" ){
              f_over_g <- dbeta(phi, h1/2, h2/2)
            } else if ( sampling_distribution == "beta" ){
              f_over_g <- dbeta(phi, h1/2, h2/2) / dbeta(phi, h1/(2 * scaling_parameter), h2/(2 * scaling_parameter))
            }
            
            # Get approximation to posterior, from which to sample
            is_distribution <- mu_pos[-(1:N_for_GP)] * sqrt(dbeta(X[-(1:N_for_GP)], h1/2, h2/2))
            
            # Sample phi values (stratified approximation)
#            inds <- sample(1:length(is_distribution), N_for_IS, replace=TRUE, prob=is_distribution)
            is_distribution_cumulative <- cumsum(is_distribution) / sum(is_distribution)
            z <- runif(N_for_IS) / N_for_IS + seq(0, 1 - 1/N_for_IS, by=1/N_for_IS)
            inds <- rep(NA, N_for_IS)
            for ( IS_ind in 1:N_for_IS ){
              inds[IS_ind] <- which.min(abs(z[IS_ind] - is_distribution_cumulative)) # re-write using apply!
            }
            phi2 <- X[-(1:N_for_GP)][inds]
            
            # Calculate whether each is in S
            alpha <- sqrt(phi2/phi_obs)
            beta <- sqrt((1 - phi2)/(1 - phi_obs))
            inS2 <- rep(NA, N_for_IS)
            for ( j in 1:N_for_IS ){
              x_phi <- x
              x_phi[(b[cp_index] - h1 + 1):(b[cp_index] + h2)] <- c(alpha[j] * x[(b[cp_index] - h1 + 1):b[cp_index]], beta[j] * x[(b[cp_index] + 1):(b[cp_index] + h2)])
              results_phi <- find_changepoints(x_phi, method=method, params=params, model=model)
              if ( length(results_phi$changepoints) >= 1 ){
                inS2[j] <- ifelse(b[cp_index] %in% results_phi$changepoints, 1, 0)
              } else {
                inS2[j] <- 0
              }
            }
            
            # Combine with initial GP samples
            phi <- c(phi, phi2)
            inS <- c(inS, inS2)
            
            # Calculate p-value estimate
            f_over_g <- c(f_over_g, dbeta(phi2, h1/2, h2/2) / (is_distribution[inds] / sum(is_distribution) * length(is_distribution)))
#            P_phi_in_S[cp_index, psi_index] <- sum(f_over_g[inS == 1])
            P_phi_in_S[cp_index, psi_index] <- sum(f_over_g[inS == 1]) / sum(f_over_g)
            p_hat[cp_index, psi_index] <- sum(f_over_g[inS == 1 & !(phi > phi_lower & phi < phi_upper)]) / sum(f_over_g[inS == 1])
#            P_both[cp_index, psi_index] <- sum(f_over_g[inS == 1 & !(phi > phi_lower & phi < phi_upper)])
#            P_phi_in_S[cp_index, psi_index] <- sum(f_over_g[inS == 1])
          }
        }
        
      } else {
        # Importance sampling, without running a GP first
        
        # Sample phi values
        if ( sampling_distribution == "uniform" ){
          phi <- runif(N_for_IS) / N_for_IS + seq(0, (N_for_IS - 1)/N_for_IS, length.out=N_for_IS)
        } else if ( sampling_distribution == "beta" ){
          z <- runif(N_for_IS) / N_for_IS + seq(0, (N_for_IS - 1)/N_for_IS, length.out=N_for_IS)
          phi <- qbeta(z, h1/(2 * scaling_parameter), h2/(2 * scaling_parameter))
        }
        
        # Calculate I(phi in S) for each phi
        alpha <- sqrt(phi/phi_obs)
        beta <- sqrt((1 - phi)/(1 - phi_obs))
        inS <- rep(NA, N_for_IS)
        for ( j in 1:N_for_IS ){
          x_phi <- x
          x_phi[(b[cp_index] - h1 + 1):(b[cp_index] + h2)] <- c(alpha[j] * x[(b[cp_index] - h1 + 1):b[cp_index]], beta[j] * x[(b[cp_index] + 1):(b[cp_index] + h2)])
          results_phi <- find_changepoints(x_phi, method=method, params=params, model=model)
          if ( length(results_phi$changepoints) >= 1 ){
            inS[j] <- ifelse(b[cp_index] %in% results_phi$changepoints, 1, 0)
          } else {
            inS[j] <- 0
          }
        }
        
        # Calculate p-value estimate
        if ( sampling_distribution == "uniform" ){
          f_over_g <- dbeta(phi, h1/2, h2/2)
        } else if ( sampling_distribution == "beta" ){
          f_over_g <- dbeta(phi, h1/2, h2/2) / dbeta(phi, h1/(2 * scaling_parameter), h2/(2 * scaling_parameter))
        }
#        P_phi_in_S[cp_index, psi_index] <- sum(f_over_g[inS == 1])
        P_phi_in_S[cp_index, psi_index] <- sum(f_over_g[inS == 1]) / sum(f_over_g)
        p_hat[cp_index, psi_index] <- sum(f_over_g[inS == 1 & !(phi > phi_lower & phi < phi_upper)]) / sum(f_over_g[inS == 1])
#        P_both[cp_index, psi_index] <- sum(f_over_g[inS == 1 & !(phi > phi_lower & phi < phi_upper)])
#        P_phi_in_S[cp_index, psi_index] <- sum(f_over_g[inS == 1])
      }
      
    }
    
  }
#  p_hat <- rowSums(P_both, na.rm=TRUE) / rowSums(P_phi_in_S, na.rm=TRUE)
  p <- rowSums(p_hat * P_phi_in_S, na.rm=TRUE) / sum(P_phi_in_S, na.rm=TRUE)
  
  return(list(changepoints=b, p=p, P_phi_in_S=P_phi_in_S, p_W=p_hat))
}
