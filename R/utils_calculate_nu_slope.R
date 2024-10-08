#create_lambda_change_in_slope <- function(n, s=1, e=n, tau=NULL, return_full=TRUE){
calculate_nu_slope <- function(n, s=1, e=n, tau=NULL, return_full=TRUE){
  
  # Calculates \nu_\tau for each possible \tau in an interval (s,e)
  
  l <- e - s + 1
  
  if ( is.null(tau) ){
    
    a <- 1/sqrt(1/6 * l * (l^2 - 1) * (1 + (e - ((s + 1):(e - 1)) + 1) * (((s + 1):(e - 1)) - s + 1) + (e - ((s + 1):(e - 1))) * (((s + 1):(e - 1)) - s)))
    b <- sqrt((e - ((s + 1):(e - 1)) + 1) * (e - ((s + 1):(e - 1))) / ((((s + 1):(e - 1)) - s) * (((s + 1):(e - 1)) - s + 1)))
    nu1 <- a * b * (3 * (((s + 1):(e - 1)) - s + 1) + e - ((s + 1):(e - 1)) - 1)
    nu2 <- a * b * (((s + 1):(e - 1)) * (e - s) + 2 * s * (((s + 1):(e - 1)) - s + 1))
    nu3 <- (-a/b) * (3 * e - 2 * ((s + 1):(e - 1)) - s + 2)
    nu4 <- (a / b) * ((s + 1):(e - 1) * (e - s) + 2 * e * (e - ((s + 1):(e - 1)) + 1))
    
    X1 <- cbind(rep(0, l), t(nu1 %*% t(s:e) - nu2))
    X1[lower.tri(X1)] <- 0
    X2 <- cbind(rep(0, l), t(nu3 %*% t(s:e) + nu4))
    X2[upper.tri(X2, diag=TRUE)] <- 0
    nu_matrix <- X1[,-1] + X2[,-1] # column i is nu for tau = i+1
    
    if ( is.null(dim(nu_matrix)) ){
      nu_matrix <- matrix(nu_matrix, ncol=1)
    }
    
    if ( return_full ){
      nu_matrix <- rbind(matrix(0, nrow=(s - 1), ncol=(e - s - 1)), nu_matrix, matrix(0, nrow=(n - e), ncol=(e - s - 1)))
    }
    
  } else {
    
    a <- 1/sqrt(1/6 * l * (l^2 - 1) * (1 + (e - tau + 1) * (tau - s + 1) + (e - tau) * (tau - s)))
    b <- sqrt((e - tau + 1) * (e - tau) / ((tau - s) * (tau - s + 1)))
    coef_lower <- a * b * (3 * (tau - s + 1) + e - tau - 1)
    cst_lower <- a * b * (tau * (e - s) + 2 * s * (tau - s + 1))
    coef_upper <- (-a/b) * (3 * e - 2 * tau - s + 2)
    cst_upper <- (a / b) * (tau * (e - s) + 2 * e * (e - tau + 1))
    
    if ( return_full ){
      nu_matrix <- matrix(0, nrow=n, ncol=1)
      nu_matrix[s:tau, 1] <- coef_lower * (s:tau) - cst_lower
      nu_matrix[(tau + 1):e, 1] <- coef_upper * ((tau + 1):e) + cst_upper
    } else {
      nu_matrix <- matrix(NA, nrow=e - s + 1, ncol=1)
      nu_matrix[1:(tau - s + 1), 1] <- coef_lower * (s:tau) - cst_lower
      nu_matrix[(tau - s + 2):(e - s + 1), 1] <- coef_upper * ((tau + 1):e) + cst_upper
    }
    
  }
  
  return(nu_matrix)
  
}




calculate_lrs_slope <- function(x, s=1, e=length(x)){
  
  l <- e - s + 1
  
  x_sums <- c(rep(0, s - 1), cumsum(x[s:e]))
  xt_sums <- c(rep(0, s - 1), cumsum((s:e) * x[s:e]))
  
  a <- b <- coef_lower <- coef_upper <- cst_lower <- cst_upper <- rep(NA, length(x))
  for ( tau in (s + 1):(e - 1) ){
    a[tau] <- 1/sqrt(1/6 * (e - s + 1) * ((e - s + 1)^2 - 1) * (1 + (e - tau + 1) * (tau - s + 1) + (e - tau) * (tau - s)))
    b[tau] <- sqrt((e - tau + 1) * (e - tau) / ((tau - s) * (tau - s + 1)))
    coef_lower[tau] <- a[tau] * b[tau] * (3 * (tau - s + 1) + (e - tau) - 1)
    cst_lower[tau] <- (-1) * a[tau] * b[tau] * (tau * (e - s) + 2 * s * (tau - s + 1))
    coef_upper[tau] <- (-1) * a[tau] / b[tau] * (3 * (e - tau) + (tau - s + 1) + 1)
    cst_upper[tau] <- a[tau] / b[tau] * (tau * (e - s) + 2 * e * (e - tau + 1))
  }
  
  nuTx <- rep(NA, length(x))
  for ( tau in (s + 1):(e - 1) ){
    nuTx[tau] <- coef_lower[tau] * xt_sums[tau] + coef_upper[tau] * (xt_sums[e] - xt_sums[tau]) + cst_lower[tau] * x_sums[tau] + cst_upper[tau] * (x_sums[e] - x_sums[tau])
  }
  
  return(nuTx)
  
}
