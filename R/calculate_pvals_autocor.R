#' Calculate p-values for change in mean model, with autocorrelation
#'
#' @param results Output of \code{find_changepoints} applied to data.
#' @param Sigma Covariance matrix of the data. This can be calculated within the function (exactly if \code{rho} and \code{sigma2} are supplied, or estimated if they are not)
#' but if the function is run multiple times where the value of \code{Sigma} is the same, it will run faster if \code{Sigma} is calculated in advance.
#' @param rho Autocorrelation parameter for AR(1) model; ignored if \code{Sigma} is supplied.
#' @param sigma2 Noise variance for AR(1) model; ignored if \code{Sigma} is supplied.
#' @param h Window size: either a positive integer (>=2) or vector of length 2 or \code{NULL}. If \code{NULL} then the changepoints either side of 
#' the changepoint of interest are used to define the window.
#' @param N Number of random \eqn{\psi}'s to generate; defaults to 1.
#' @param est_sigma_pars List of parameters to be supplied to the function \code{est_rho}, in order to estimate the variance and autocorrelation parameter.
#' This will be ignored if either \code{Sigma} or both \code{rho} and \code{sigma2} are supplied.
#' @param cp_bound If \code{TRUE}, then if there is an estimated changepoint within the window that is fixed under the null hypothesis,
#' this changepoint will be used as the boundary of the window.
#' @param gamma Number in \code{(0, 1]}. If \code{h = NULL}, then \eqn{h} is taken to be \code{gamma} times the distance between the changepoint
#' of interest and the estimated changepoints on either side. Ignored if \code{h} specified; defaults to \code{1}.
#' @param num_pvals Number of p-values to compute. If \code{NULL}, the function will calculate p-values for all detected changepoints.
#' 
#' @description Calculate p-values for detected changepoints, for change in mean or slope, where there is autocorrelated noise. Changepoints can be detected using
#' the CUSUM statistic with either binary segmentation, wild binary segmentation, or narrowest over threshold.
#'
#' @details 
#' The model is \eqn{X_t = \mu_t + Z_t}, where \eqn{\mu_t} is constant except at changepoints and \eqn{Z_t \sim N(0, \Sigma)}.
#' Either the covariance matrix \eqn{\Sigma} can be supplied to the function or, if \eqn{Z_t} follows an AR(1) model of the form
#' \eqn{Z_t = \rho Z_{t-1} + \sigma \epsilon_t} with eqn{\epsilon_t \sim N(0, 1)}, then \eqn{\rho} and \eqn{\sigma^2} can be supplied instead.
#' If neither of \code{Sigma} or \code{rho} is specified, then the function will estimate \eqn{rho} and \eqn{\sigma^2} under the assumption
#' that the data follows an AR(1) model.
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
#' @return A list:
#' \itemize{
#' \item \code{p_value} A vector containing p-values for each estimated changepoint.
#' \item \code{P_both} A matrix containing values of \eqn{Pr(\phi \in S \& |\phi| > |\phi_{obs}|)} for each estimated changepoint and sampled \eqn{\psi} value.
#' \item \code{P_phi_in_S} A matrix containing values of \eqn{Pr(\phi \in S)} for each estimated changepoint and sampled \eqn{\psi} value.
#' }
#'
#' @export
#'
#' @examples 
#' set.seed(1)
#' x <- arp(100, 0.6, c(rep(0, 50), rep(2, 50)))
#' results <- find_changepoints(x, method="bs", model="mean")
#' calculate_pvals_autocor(results, rho=0.6, sigma2=1, h=20)
#' calculate_pvals_autocor(results, h=20) # using estimated variances
#' 
calculate_pvals_autocor <- function(results, Sigma=NULL, rho=NULL, sigma2=NULL, h=NULL, N=1, est_sigma_pars=list(K=10, rho_grid_size=0.01), cp_bound=TRUE, gamma=1, num_pvals=NULL){
  
  x <- results$x
  n <- length(x)
  
  b <- results$changepoints
  d <- results$results$d
  method <- results$params$method
  model <- results$params$model
  
  if ( !(model %in% c("mean", "slope")) ){
    stop("Incorrect model.")
  }
  
  if ( length(b) == 0 ){
    stop("No changepoints detected.")
  } else if ( is.na(b[1]) ){
    stop("No changepoints detected.")
  } else if ( is.null(num_pvals) ){
    num_pvals <- length(b)
  } else if ( length(b) < num_pvals ){
    print(paste0("Only ", length(b), " changepoints were detected."))
    num_pvals <- length(b)
  }
  
  p_value <- rep(NA, num_pvals)
  P_both <- P_phi_in_S <- matrix(NA, nrow=num_pvals, ncol=N)
  
  if ( is.null(Sigma) ){
    if ( is.null(rho) || is.null(sigma2) ){
      if ( model == "mean" ){
        sig_ests <- do.call(est_rho, c(est_sigma_pars, x=list(x)))
        rho <- sig_ests$rho
        sigma2 <- sig_ests$sigma2
      } else if ( model == "slope" ){
        sig_ests <- do.call(est_rho, c(est_sigma_pars, x=list(x), model="slope"))
        rho <- sig_ests$rho
        sigma2 <- sig_ests$sigma2
      }
    }
    Sigma <- diag((1 - rho^seq(2, 2*n, by=2))/(1 - rho^2))
    for ( i in 1:(n - 1) ){
      Sigma[(i + 1):n, i] <- Sigma[i, (i + 1):n] <- rho^seq(1:(n - i)) * Sigma[i, i]
    }
    Sigma <- sigma2 * Sigma
  }
  if ( N > 1 ){
    svd_Sigma <- svd(Sigma)
    D <- svd_Sigma$d
    V <- svd_Sigma$v
  }

  for ( cp_index in 1:num_pvals ){
    
    # Calculate h
    hvals <- get_h1_h2(h, b, n, b[cp_index], cp_bound, gamma)
    h1 <- hvals[1]
    h2 <- hvals[2]
    nh <- h1 + h2
    
    # Calculate nu
    if ( model == "mean" ){
      nu <- c(rep(0, b[cp_index] - h1), rep(1/h1, h1), rep(-1/h2, h2), rep(0, n - b[cp_index] - h2))
    } else if ( model == "slope" ){
      nu <- calculate_nu_slope(n, b[cp_index] - h1 + 1, b[cp_index] + h2, b[cp_index], return_full=TRUE)
      nu <- nu / sqrt(sum(nu^2))
    }
    phi_obs <- c(t(nu[(b[cp_index] - h1 + 1):(b[cp_index] + h2)]) %*% x[(b[cp_index] - h1 + 1):(b[cp_index] + h2)])
    
    # Calculate variance of phi
    phi_var <- c(t(nu[(b[cp_index] - h1 + 1):(b[cp_index] + h2)]) %*% 
                        Sigma[(b[cp_index] - h1 + 1):(b[cp_index] + h2), (b[cp_index] - h1 + 1):(b[cp_index] + h2)] %*% 
                        nu[(b[cp_index] - h1 + 1):(b[cp_index] + h2)])
    phi_sd <- sqrt(phi_var)

    # Calculate S_all
    if ( length(b) >= 1 & !is.null(h) ){
      S_all <- calculate_S_autocor(x, results=results, Sigma=Sigma, nu=nu, phi_obs=phi_obs, phi_var=phi_var, first_cp_only=TRUE, cp_of_interest=cp_index)
    } else {
      S_all <- calculate_S_autocor(x, results=results, Sigma=Sigma, nu=nu, phi_obs=phi_obs, phi_var=phi_var)
    }
    
    # Calculate S
    S <- get_S(S_all, h, b, tau=b[cp_index])
    
    # Calculate P(phi \in S)
    P_phi_in_S[cp_index, 1] <- sum(pnorm(S[,2] / phi_sd) - pnorm(S[,1] / phi_sd))
    
    # Calculate P(phi > |phi_obs| & phi \in S)
    P_both[cp_index, 1] <- 0
    if ( nrow(S) >= 1 ){
      for ( i in 1:nrow(S) ){
        if ( S$upper_lim[i] > abs(phi_obs) ){
          P_both[cp_index, 1] <- P_both[cp_index, 1] + pnorm(S$upper_lim[i] / phi_sd) - pnorm(max(c(abs(phi_obs), S$lower_lim[i])) / phi_sd)
        }
        if ( S$lower_lim[i] < (-1) * abs(phi_obs) ){
          P_both[cp_index, 1] <- P_both[cp_index, 1] + pnorm(min(c(-abs(phi_obs), S$upper_lim[i])) / phi_sd) - pnorm(S$lower_lim[i] / phi_sd)
        }
      }
    }
      
    if ( N > 1 ){
      nu_tilde <- diag(sqrt(D)) %*% t(V) %*% nu
      x_tilde <- diag(1/sqrt(D)) %*% t(V) %*% x
      
      W <- rbind(cbind(matrix(0, 1, b[cp_index] - h1), matrix(1, 1, h1 + h2), matrix(0, 1, n - b[cp_index] - h2)),
                 cbind(diag(b[cp_index] - h1), matrix(0, b[cp_index] - h1, n - b[cp_index] + h1)),
                 cbind(matrix(0, n - b[cp_index] - h2, b[cp_index] + h2), diag(n - b[cp_index] - h2)))
      W_tilde <- W %*% V %*% diag(1/sqrt(D))
      V_W <- svd(W_tilde)$v
      
      Z_tilde <- diag(n) - crossprod(t(V_W)) - 1/phi_var * crossprod(t(nu_tilde))
      U_tilde <- try(svd::propack.svd(Z_tilde, neig=h1 + h2 - 2)$u)
      if ( class(U_tilde)[1] == "try-error" ){
        U_tilde <- try(svd::propack.svd(Z_tilde)$u[,1:(h1 + h2 - 2)])
      }
      if ( class(U_tilde)[1] == "try-error" ){
        U_tilde <- try(svd(Z_tilde)$u[,1:(h1 + h2 - 2)])
      } else if ( sum(is.na(U_tilde)) > 0 ){
        U_tilde <- try(svd(Z_tilde)$u[,1:(h1 + h2 - 2)])
      }
      
      if ( class(U_tilde)[1] != "try-error" ){
        if( dim(U_tilde)[2] == h1 + h2 - 2 ){
          if ( sum(is.na(U_tilde)) == 0 ){
      
        psi_mat <- matrix(rnorm((N - 1) * (h1 + h2 - 2)), nrow=N - 1)
        
        x0 <- V %*% diag(sqrt(D)) %*% crossprod(t(V_W)) %*% diag(1/sqrt(D)) %*% t(V) %*% x + 1/c((t(nu) %*% Sigma %*% nu)) * Sigma %*% nu * phi_obs
        psi_coef <- V %*% diag(sqrt(D)) %*% U_tilde
        
        for ( psi_index in 2:N ){
          psi <- psi_mat[psi_index - 1,]
          x_psi <- x0 + psi_coef %*% psi
          results_psi <- find_changepoints(x_psi, method=method, params=within(results$params, rm(method, model)))
          
          S_all <- calculate_S_autocor(x_psi, nu=nu, Sigma=Sigma, results=results_psi, phi_obs=phi_obs, phi_var=phi_var)
          S <- get_S(S_all, h, b, tau=b[cp_index])
          
          P_phi_in_S[cp_index, psi_index] <- sum(pnorm(S[,2] / phi_sd) - pnorm(S[,1] / phi_sd))
          P_both[cp_index, psi_index] <- 0
          if ( nrow(S) >= 1 ){
            for ( k in 1:nrow(S) ){
              if ( S$upper_lim[k] > abs(phi_obs) ){
                P_both[cp_index, psi_index] <- P_both[cp_index, psi_index] +
                  pnorm(S$upper_lim[k] / phi_sd) - pnorm(max(c(abs(phi_obs), S$lower_lim[k])) / phi_sd)
              }
              if ( S$lower_lim[k] < (-1) * abs(phi_obs) ){
                P_both[cp_index, psi_index] <- P_both[cp_index, psi_index] +
                  pnorm(min(c(-abs(phi_obs), S$upper_lim[k])) / phi_sd) - pnorm(S$lower_lim[k] / phi_sd)
              }
            }
          }
        }

      }
      
    }}
    }
    
    # Calculate p-value
    p_value[cp_index] <- sum(P_both[cp_index,]) / sum(P_phi_in_S[cp_index,])
    
    if ( P_phi_in_S[cp_index, 1] < 10^(-12) ){
      warning("The probability of the conditioning event is very small, p-value estimate may be unstable.")
    }
    
  }
  
  return(list(p_value=p_value, P_both=P_both, P_phi_in_S=P_phi_in_S))
  
}