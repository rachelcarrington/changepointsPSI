# Calculate vector of CUSUM statistics for data in terms of phi.
# Used inside \code{calculate_pvals}. This file contains \code{cusum_phi_vec} and \code{cusum_phi_vec_var}.

cusum_phi_vec <- function(x, nu, nu2=NULL, phi_obs=NULL, s=1, e=length(x), autocor=FALSE, icss=FALSE, C0=NULL, C0_nu=NULL, Sigma=NULL){
  
  # Returns \code{length(x)} x 2 matrix such that cusum_phi_vec(x, nu)^T %*% c(1, phi) = cusum(x_phi(x, nu, phi)).
  
  # x Vector of data
  # nu Contrast vector: the test statistic is nu^T x
  # nu2 Value of ||nu||_2^2
  # phi_obs Observed value of test statistic (nu^T x)
  # s Starting point for calculating CUSUM statistics
  # e Ending point for calculating CUSUM statistics
  # autocor Whether there is autocorrelation in the data
  # icss Whether to use iterated cumulative sum of squares rather than CUSUM
  # C0 Vector of cumulative sums of x
  # C0_nu Vector of cumulative sums of nu
  # Sigma Covariance matrix of x (if autocor = TRUE)

  if ( is.null(nu2) ){
    if ( !autocor ){
      nu2 <- sum(nu^2)
    } else {
      nu2 <- c(t(nu) %*% Sigma %*% nu)
    }
  }

  if ( is.null(phi_obs) ){
    phi_obs <- c(t(nu) %*% x)
  }
  
  if ( !autocor ){
    cusum_nu <- cusum(nu, s, e, cumsums=C0_nu, icss=icss)
  } else {
    cusum_nu <- cusum(Sigma %*% nu, s, e, cumsums=C0_nu, icss=icss)
  }
  cst <- cusum(x, s, e, cumsums=C0, icss=icss) - phi_obs/nu2 * cusum_nu
  coef <- cusum_nu / nu2

  return(cbind(cst, coef))

}


cusum_phi_vec_var <- function(x, tau, h1, h2, s=1, e=NULL, return_full=FALSE, cp_bound=TRUE, gamma=1, icss=FALSE, C0=NULL){

  n <- length(x)

  if ( is.null(e) ){
    e <- n
  }
  
  if ( is.null(C0) ){
    C0 <- cumsum(x^2)
  }

  if ( tau - h1 == 0 ){
    S0 <- C0[tau + h2]
    S1 <- C0[tau]
  } else {
    S0 <- C0[tau + h2] - C0[tau - h1]
    S1 <- C0[tau] - C0[tau - h1]
  }
  S2 <- C0[tau + h2] - C0[tau]
  
  a0 <- c(rep(1, tau - h1), rep(0, h1), rep(S0/S2, h2), rep(1, n - tau - h2))
  a1 <- c(rep(0, tau - h1), rep(S0/S1, h1), rep(-S0/S2, h2), rep(0, n - tau - h2))
  
  cst <- coef <- rep(NA, n - 1)
  for ( j in s:(e - 1) ){
    cst[j] <- sqrt((j - s + 1) * (e - j)/(e - s + 1)) * (mean(a0[s:j] * x[s:j]^2) - mean(a0[(j + 1):e] * x[(j + 1):e]^2))
    coef[j] <- sqrt((j - s + 1) * (e - j)/(e - s + 1)) * (mean(a1[s:j] * x[s:j]^2) - mean(a1[(j + 1):e] * x[(j + 1):e]^2))
  }
  cst[abs(cst) < 10^(-10)] <- 0
  coef[abs(coef) < 10^(-10)] <- 0
  
  if ( !(return_full) ){
    cst <- cst[s:(e - 1)]
    coef <- coef[s:(e - 1)]
  }
  
  return(cbind(cst, coef))

}