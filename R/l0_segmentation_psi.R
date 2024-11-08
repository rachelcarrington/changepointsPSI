#' Post-selection inference for the change in mean model, using L0 segmentation to estimate changepoints.
#'
#' @description Post-selection inference for L0 segmentation. This uses the functions \code{changepoint_estimates} and
#' \code{changepoint_inference} from the package \code{ChangepointInference}.
#'
#' @param x Numeric vector of data.
#' @param lambda Changepoint detection threshold.
#' @param h Window size.
#' @param N Number of \eqn{\psi} samples to take (defaults to 1).
#' @param sigma2 Variance of \code{x}.
#' @param sig Tuning parameter.
#' @param include_original Logical; whether to include observed value as \eqn{\psi} in place as one of the random samples;
#' defaults to \code{TRUE}.
#' @param num_pvals Maximum number of p-values to calculate.
#' @param cp_bound If \code{TRUE}, then if there is an estimated changepoint within the window that is fixed under the null hypothesis,
#' this changepoint will be used as the boundary of the window.
#' @param gamma Value in \code{(0, 1]}. If \code{h = NULL}, then \eqn{h} is taken to be \code{gamma} times the distance between the changepoint
#' of interest and the estimated changepoints on either side. Ignored if \code{h} specified; defaults to \code{1}.
#'
#' @return A list:
#' \itemize{
#' \item \code{b} Vector of changepoints
#' \item \code{p_value} Vector of p-values
#' \item \code{p_value_orig} Vector of p-values obtained using fixed \eqn{\psi = \psi_{obs}}
#' \item \code{phi_obs} The observed value of the test statistic.
#' \item \code{P_both} Matrix containing values of \eqn{Pr(|\phi| > |\phi_{obs}| \& \phi \in S)}
#' \item \code{P_phi_in_S} Matrix containing values of \eqn{Pr(\phi \in S)}
#' \item \code{P_both_orig} Vector containing values of \eqn{Pr(|\phi| > |\phi_{obs}| \& \phi \in S)} for fixed \eqn{\psi = \psi_{obs}}
#' \item \code{P_phi_in_S_orig} Vector containing values of \eqn{Pr(\phi \in S)} for fixed \eqn{\psi = \psi_{obs}}
#' }
#'
#' @importFrom ChangepointInference changepoint_estimates
#' @importFrom ChangepointInference changepoint_inference
#' @export
#'
#' @examples
#' set.seed(100)
#' x <- rnorm(100) + c(rep(1,40), rep(-1,20), rep(1,40))
#' l0_segmentation_psi(x, 4, 10, 10)
#'
l0_segmentation_psi <- function(x, lambda, h, N=1, sigma2=NULL, sig=4, include_original=TRUE, num_pvals=NULL, cp_bound=TRUE, gamma=1){
  
  # Requires package ChangepointInference
  
  n <- length(x)
  
  if ( is.null(sigma2) ){
    sigma2 <- mad(x)
  }
  
  fit <- ChangepointInference::changepoint_estimates(x, "L0", lambda)
  b <- fit$change_pts
  if ( is.null(num_pvals) ){
    num_pvals <- length(fit$change_pts)
  } else if ( num_pvals > length(fit$change_pts) ){
    num_pvals <- length(fit$change_pts)
  }
  
  if ( length(b) >= 1 ){
    
    if ( include_original ){
      
      # Original p-values (N = 1)
      P_both_orig <- P_phi_in_S_orig <- rep(0, num_pvals)
      pvals_orig <- rep(NA, num_pvals)
      l0_results <- ChangepointInference::changepoint_inference(x, "L0-fixed", tuning_parameter=lambda, window_size=h, sig=sig, return_conditioning_set=TRUE)
      
      for ( cp_index in 1:num_pvals ){
        hvals <- get_h1_h2(h, b, n, b[cp_index], cp_bound, gamma)
        h1 <- hvals[1]
        h2 <- hvals[2]
        nu <- c(rep(0, b[cp_index] - h1), rep(1/h1, h1), rep(-1/h2, h2), rep(0, n - b[cp_index] - h2))
        nu2 <- sum(nu^2)
        phi_sd <- sqrt(nu2 * sigma2)
        phi_obs <- c(t(nu) %*% x)
        
        S <- l0_results$conditioning_sets[[cp_index]]
        S2 <- S[S$contained == 1, ]
        
        # Calculate P(phi \in S)
        P_phi_in_S_orig[cp_index] <- sum(pnorm(S2$max_mean / phi_sd) - pnorm(S2$min_mean / phi_sd))
        
        # Calculate P(|phi| > |phi_obs| & phi \in S)
        for ( i in 1:nrow(S2) ){
          if ( S2$max_mean[i] > abs(phi_obs) ){
            P_both_orig[cp_index] <- P_both_orig[cp_index] + pnorm(S2$max_mean[i] / phi_sd) - pnorm(max(c(abs(phi_obs), S2$min_mean[i])) / phi_sd)
          }
          if ( S2$min_mean[i] < (-1) * abs(phi_obs) ){
            P_both_orig[cp_index] <- P_both_orig[cp_index] + pnorm(min(c(-abs(phi_obs), S2$max_mean[i])) / phi_sd) - pnorm(S2$min_mean[i] / phi_sd)
          }
        }
        pvals_orig[cp_index] <- P_both_orig[cp_index] / P_phi_in_S_orig[cp_index]
      }
      
      N2 <- N - 1
      
    } else {
      
      pvals_orig <- P_both_orig <- P_phi_in_S_orig <- NA
      N2 <- N
      
    }
    
    # New p-values
    
    if ( N2 >= 1 ){
      
      alpha <- rep(1, h1 + h2)
      alpha2 <- sum(alpha^2)
      Z <- diag(h1 + h2) - 1/alpha2 * crossprod(t(alpha)) - 1/(h1 + h2) * crossprod(t(c(rep(1/h1, h1), rep(-1/h2, h2))))
      U <- svd(Z)$u[,1:(h1 + h2 - 2)]
      
      if ( is.null(num_pvals) ){
        num_pvals <- length(b)
      } else {
        num_pvals <- min(c(num_pvals, length(b)))
      }
      
      p_val <- rep(NA, num_pvals)
      P_both <- P_phi_in_S <- matrix(NA, nrow=num_pvals, ncol=N2)
      
      for ( cp_index in 1:num_pvals ){
        
#        if ( b[cp_index] >= h & b[cp_index] <= n - h - 1 ){
        
        nu <- c(rep(0, b[cp_index] - h1), rep(1/h1, h1), rep(-1/h2, h2), rep(0, n - b[cp_index] - h2))
        nu2 <- sum(nu^2)
        phi_sd <- sqrt(nu2 * sigma2)
        phi_obs <- c(t(nu) %*% x)
        
        Psi <- matrix(rnorm(N2 * (h1 + h2 - 2), sd=sqrt(sigma2)), nrow=N2)
        for ( psi_index in 1:N2 ){
          x_psi <- x
          x_psi[(b[cp_index] - h1 + 1):(b[cp_index] + h2)] <- x_psi[(b[cp_index] - h1 + 1):(b[cp_index] + h2)] + U %*% Psi[psi_index,]
          
          fit_new <- ChangepointInference::changepoint_estimates(x_psi, "L0", lambda)
          b_new <- fit_new$change_pts
          
          # I think this works (we need to make sure that b[cp_index] %in% b_new)
          phi <- 10
          while ( !(b[cp_index] %in% b_new) ){
            x_psi <- calculate_x_phi(x_psi, nu, phi)
            phi <- phi + 10
            b_new <- ChangepointInference::changepoint_estimates(x_psi, "L0", lambda)$change_pts
          }
          
          results_l0 <- ChangepointInference::changepoint_inference(x_psi, 'L0-fixed', lambda, window_size=h, sig=sig, return_conditioning_sets=TRUE)
          
          # Make sure we select the right changepoint!
          S <- results_l0$conditioning_sets[[(1:length(results_l0$change_pts))[results_l0$change_pts == b[cp_index]]]]
          
          # Keep rows of S for which b[cp_index] %in% b ("contained" tells us whether b_new[1] %in% b)
          S2 <- S[S$contained == 1,]
          
          # Calculate P(phi \in S)
          P_phi_in_S[cp_index, psi_index] <- sum(pnorm(S2$max_mean / phi_sd) - pnorm(S2$min_mean / phi_sd))
          
          # Calculate P(phi > |phi_obs| & phi \in S)
          P_both[cp_index, psi_index] <- 0
          for ( i in 1:nrow(S2) ){
            if ( S2$max_mean[i] > abs(phi_obs) ){
              P_both[cp_index, psi_index] <- P_both[cp_index, psi_index] + pnorm(S2$max_mean[i] / phi_sd) - pnorm(max(c(abs(phi_obs), S2$min_mean[i])) / phi_sd)
            }
            if ( S2$min_mean[i] < (-1) * abs(phi_obs) ){
              P_both[cp_index, psi_index] <- P_both[cp_index, psi_index] + pnorm(min(c(-abs(phi_obs), S2$max_mean[i])) / phi_sd) - pnorm(S2$min_mean[i] / phi_sd)
            }
          }
          
        } # end loop over psi
        
        # Calculate p-value estimates
        if ( include_original ){
          p_val[cp_index] <- (sum(P_both[cp_index,]) + P_both_orig[cp_index]) / (sum(P_phi_in_S[cp_index,]) + P_phi_in_S_orig[cp_index])
        } else {
          p_val[cp_index] <- sum(P_both[cp_index,]) / sum(P_phi_in_S[cp_index,])
        } 
      
      } # end loop over b
      
    } else {
      
      P_both <- P_both_orig
      P_phi_in_S <- P_phi_in_S_orig
      p_val <- pvals_orig
      
    }
    
    
  } else {
    
    p_val <- pvals_orig <- P_both <- P_phi_in_S <- P_phi_in_S_orig <- P_both_orig <- phi_obs <- NA
    print("No changepoints detected.")
    
  }
  
  return(list(b=b, p_value=p_val, p_value_orig=pvals_orig, phi_obs=phi_obs, P_both=P_both, P_phi_in_S=P_phi_in_S, P_both_orig=P_both_orig,
                P_phi_in_S_orig=P_phi_in_S_orig))
}

