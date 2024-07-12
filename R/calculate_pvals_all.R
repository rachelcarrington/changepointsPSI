#' Calculate p-values for changepoints
#'
#' @description Calculate p-values for detected changepoints, for change in mean, slope, or variance. Changepoints can be detected using
#' the CUSUM statistic with either binary segmentation, wild binary segmentation, or narrowest over threshold.
#'
#' @param results Output of \code{find_changepoints}.
#' @param h Window size: the null hypothesis is that there are no changepoints within a window \code{h} of the changepoint of interest. This can be a positive integer (>=2) 
#' or a vector of length 2 (if the window is asymmetric) or \code{NULL} (if the window is defined by the changepoints either side of the changepoint of interest).
#' @param N Number of random \eqn{\psi}'s to generate; defaults to 1.
#' @param num_pvals Maximum number of p-values to calculate; if set to NA, it will default to \code{length(x) - 1}. P-values are calculated for changepoints
#' in the order of detection for iterative algorithms, otherwise they are calculated in the numeric order of changepoint locations.
#' @param sigma2 Variance of \code{x}, for piecewise constant/linear mean models. If not given, it will be estimated using median absolute deviation.
#' @param include_original Whether to include the \eqn{\psi} value corresponding to the observed data in place of one of the random samples. Defaults to \code{TRUE}.
#' @param cp_bound If \code{TRUE}, then if there is an estimated changepoint within the window that is fixed under the null hypothesis,
#' this changepoint will be used as the boundary of the window.
#' @param gamma Value in \code{(0, 1]}. If \code{h = NULL}, then \eqn{h} is taken to be \code{gamma} times the distance between the changepoint
#' of interest and the estimated changepoints on either side. Ignored if \code{h} specified; defaults to \code{1}.
#' @param return_probs Whether to return the values of \eqn{Pr(\phi \in S \& |\phi| > |\phi_{obs}|)} and \eqn{Pr(\phi \in S)}
#' for each \eqn{\psi} as part of the output. Defaults to \code{FALSE}.
#'
#' @details
#' Current options for changepoint models are change in mean (\code{model = "mean"}), change in slope (\code{model = "slope"}), and change in variance
#' (\code{model = "var"}). Either binary segmentation (\code{"bs"}), wild binary segmentation (\code{"wbs"}) or narrowest-over-threshold (\code{"not"}) can be used to estimate changepoints.
#' For the change in variance model, the CUSUM statistic must be used as the loss function (\code{loss = "cusum"}): otherwise, the p-values cannot be calculated exactly.
#' In this case p-values should be estimated using the function \code{calculate_pvals_var}.
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
#' @return A list.
#' \itemize{
#' \item \code{p_value} A vector of p-values
#' \item \code{P_both} If \code{return_probs=TRUE}, a matrix containing values of \eqn{Pr(\phi \in S \& |\phi| > |\phi_{obs}|)} for each
#' \eqn{\psi}; otherwise \code{NA}.
#' \item \code{P_phi_in_S} If \code{return_probs=TRUE}, a matrix containing values of \eqn{Pr(\phi \in S)} for each \eqn{\psi};
#' otherwise \code{NA}.
#' }
#'
#' @export
#'
#' @examples
#' set.seed(100)
#' x <- c(rep(1,50), rep(-1,50), rep(1,50), rep(-1,50)) + rnorm(200)
#' results <- find_changepoints(x, method="bs", params=list(threshold=4))
#' calculate_pvals_all(results, h=20, return_probs=TRUE)
#' calculate_pvals_all(results, h=20, N=5, return_probs=TRUE)
#'
#' x <- c(rep(0, 100), seq(0, 10, length.out=100)) + rnorm(200)
#' results <- find_changepoints(x, method="bs", model="slope")
#' calculate_pvals_all(results, h=20, sigma2=1, return_probs=TRUE)
#' calculate_pvals_all(results, h=20, N=5, sigma2=1, return_probs=TRUE)
#'
#' x <- c(rnorm(100), rnorm(100, sd=2))
#' results <- find_changepoints(x, method="bs", params=list(threshold=10, loss="cusum"), model="var")
#' calculate_pvals_all(results, h=20, return_probs=TRUE)
#' calculate_pvals_all(results, h=20, N=5, return_probs=TRUE)
#'
calculate_pvals_all <- function(results, h=NULL, N=1, num_pvals=NULL, sigma2=NULL, include_original=TRUE, cp_bound=TRUE, gamma=1, return_probs=FALSE){

  model <- results$params$model
  exact_results <- TRUE
  if ( model == "var" ){
    if ( results$params$loss == "lrs" || results$params$method == "pelt" ){
      warning("P-values cannot be calculated exactly, calculating approximate p-values.")
      exact_results <- FALSE
      return(calculate_pvals_var(results, h, num_pvals=num_pvals, NW=N))
    }
  }
  
  if ( exact_results ){
    x <- results$x
    model <- results$params$model
    method <- results$params$method
    cp_params <- within(results$params, rm(model, method))
    b <- results$changepoints
    d <- results$results$d
    n <- length(x)
    
    # Check how many changepoints are detected, and set number of p-values to calculate
    if ( length(b) == 0 ){
      stop("No changepoints detected")
    } else if ( is.null(num_pvals) ){
      num_pvals <- length(b)
    } else if ( length(b) < num_pvals ){
      num_pvals <- length(b)
    }
    
    # Estimate variance of the data, if not given
    if ( is.null(sigma2) & !(model == "var") ){
      sigma2 <- mad(x)
    }
    
    # Create vector of p-values
    p_value <- rep(NA, num_pvals)
    P_both <- P_phi_in_S <- matrix(0, nrow=num_pvals, ncol=N)
    
    for ( cp_index in 1:num_pvals ){
      
      # Calculate h
      hvals <- get_h1_h2(h, b, n, b[cp_index], cp_bound, gamma)
      h1 <- hvals[1]
      h2 <- hvals[2]
      nh <- h1 + h2

      # Calculate Psi
      if ( model == "mean" ){
        nu <- c(rep(0, b[cp_index] - h1), rep(1/h1, h1), rep(-1/h2, h2), rep(0, n - b[cp_index] - h2))
        nu2 <- 1/h1 + 1/h2
        phi_obs <- c(t(nu) %*% x) # observed value of test statistic
        phi_sd <- sqrt(sigma2 * nu2) # variance of phi
        psi_vals <- get_psi_vals(x, h1, h2, N, b[cp_index], include_original, model, sigma2=sigma2, nu=nu)
        Psi <- psi_vals$Psi
        X <- psi_vals$X
        U <- psi_vals$U
      } else if ( model == "var" ){
        C0 <- sum(x[(b[cp_index] - h1 + 1):(b[cp_index] + h2)]^2)
        C1 <- sum(x[(b[cp_index] - h1 + 1):b[cp_index]]^2)
        phi_obs <- C1 / C0
        phi_star <- qbeta(1 - pbeta(phi_obs, h1/2, h2/2), h1/2, h2/2)
        phi_upper <- max(phi_obs, phi_star)
        phi_lower <- min(phi_obs, phi_star)
        Psi <- get_psi_vals(x, h1, h2, N, b[cp_index], include_original, model, C0=C0, C1=C1)$Psi
      } else if ( model == "slope" ){
        nu <- calculate_nu_slope(n, b[cp_index] - h1 + 1, b[cp_index] + h2, b[cp_index], return_full=TRUE)
        nu <- nu / sqrt(sum(nu^2))
        phi_obs <- c(t(nu) %*% x)
        phi_sd <- sqrt(sigma2) # variance of phi
        psi_vals <- get_psi_vals(x, h1, h2, N, b[cp_index], include_original, model, sigma2=sigma2, nu=nu)
        Psi <- psi_vals$Psi
        X <- psi_vals$X
        U <- psi_vals$U
      }
      
      # Loop over psi values
      for ( iter in 1:N ){
        
        # Calculate x'(psi)
        if ( model == "mean" ){
          x_psi <- X[,1] + nu * phi_obs / nu2 # this is necessary b/c when calculating S it assumes we need to subtract it
          x_psi[(b[cp_index] - h1 + 1):(b[cp_index] + h2)] <- x_psi[(b[cp_index] - h1 + 1):(b[cp_index] + h2)] + U %*% t(Psi[iter,,drop=FALSE])
        } else if ( model == "var" ){
          x_psi <- x
          x_psi[(b[cp_index] - h1 + 1):(b[cp_index] + h2)] <- Psi[iter,]
        } else if ( model == "slope" ){
          x_psi <- X[,1] + nu * phi_obs # this is necessary b/c when calculating S it assumes we need to subtract it
          x_psi[(b[cp_index] - h1 + 1):(b[cp_index] + h2)] <- x_psi[(b[cp_index] - h1 + 1):(b[cp_index] + h2)] + U %*% t(Psi[iter,,drop=FALSE])
        }
        
        # Find changepoints
        results_psi <- find_changepoints(x_psi, method, cp_params, model)
        
        # Calculate S
        if ( model == "mean" || model == "slope" ){
          S_all <- calculate_S(x_psi, results_psi, nu=nu, phi_obs=phi_obs)
        } else if ( model == "var" ){
          S_all <- calculate_S(x_psi, results_psi, tau=b[cp_index], h=c(h1, h2), phi_obs=phi_obs)
        }
        S <- get_S(S_all, h, b, tau=b[cp_index])
        
        # Calculate P(phi in S) and P(phi in S and C)
        if ( model == "mean" ){
          P_phi_in_S[cp_index, iter] <- sum(pnorm(S[,2] / phi_sd) - pnorm(S[,1] / phi_sd))
          
          P_both[cp_index, iter] <- 0
          if ( nrow(S) >= 1 ){
            for ( i in 1:nrow(S) ){
              if ( S$upper_lim[i] > abs(phi_obs) ){
                P_both[cp_index, iter] <- P_both[cp_index, iter] + pnorm(S$upper_lim[i] / phi_sd) - pnorm(max(c(abs(phi_obs), S$lower_lim[i])) / phi_sd)
              }
              if ( S$lower_lim[i] < (-1) * abs(phi_obs) ){
                P_both[cp_index, iter] <- P_both[cp_index, iter] + pnorm(min(c(-abs(phi_obs), S$upper_lim[i])) / phi_sd) - pnorm(S$lower_lim[i] / phi_sd)
              }
            }
          }
          
        } else if ( model == "slope" ){
          P_phi_in_S[cp_index, iter] <- sum(pnorm(S[,2] / phi_sd) - pnorm(S[,1] / phi_sd))
          
          P_both[cp_index, iter] <- 0
          if ( nrow(S) >= 1 ){
            for ( i in 1:nrow(S) ){
              if ( S$upper_lim[i] > abs(phi_obs) ){
                P_both[cp_index, iter] <- P_both[cp_index, iter] + pnorm(S$upper_lim[i] / phi_sd) - pnorm(max(c(abs(phi_obs), S$lower_lim[i])) / phi_sd)
              }
              if ( S$lower_lim[i] < (-1) * abs(phi_obs) ){
                P_both[cp_index, iter] <- P_both[cp_index, iter] + pnorm(min(c(-abs(phi_obs), S$upper_lim[i])) / phi_sd) - pnorm(S$lower_lim[i] / phi_sd)
              }
            }
          }
          
        } else if ( model == "var" ){
          P_phi_in_S[cp_index, iter] <- sum(pbeta(S$upper_lim, h1/2, h2/2) - pbeta(S$lower_lim, h1/2, h2/2))
          
          P_both[cp_index, iter] <- 0
          if ( P_phi_in_S[cp_index, iter] != 0 ){
            for ( i in 1:nrow(S) ){
              if ( phi_lower > S$lower_lim[i] ){
                P_both[cp_index, iter] <- P_both[cp_index, iter] + 
                  pbeta(min(S$upper_lim[i], phi_lower), h1/2, h2/2) - pbeta(S$lower_lim[i], h1/2, h2/2)
              }
              if ( phi_upper < S$upper_lim[i] ){
                P_both[cp_index, iter] <- P_both[cp_index, iter] + 
                  pbeta(S$upper_lim[i], h1/2, h2/2) - pbeta(max(S$lower_lim[i], phi_upper), h1/2, h2/2)
              }
            }
          }
          
        } 
      } # end loop over Psi
      
      p_value[cp_index] <- sum(P_both[cp_index,]) / sum(P_phi_in_S[cp_index,])
      
      if ( P_phi_in_S[cp_index, 1] < 10^(-12) ){
        warning("The probability of the conditioning event is very small, p-value estimate may be unstable.")
      }
      
    } # end loop over b
    
    if ( !return_probs ){
      P_both <- P_phi_in_S <- NA
    }
    
    return(list(p_value=p_value, P_both=P_both, P_phi_in_S=P_phi_in_S))
    
  }
  
}
