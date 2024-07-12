# Generate Psi within calculate_pvals.

get_psi_vals <- function(x, h1, h2, N, tau, include_original, model, sigma2=1, nu=NULL, C0=NULL, C1=NULL){

  n <- length(x)
  nh <- h1 + h2

  if ( model == "mean" ){

    nu <- c(rep(1/h1, h1), rep(-1/h2, h2))
    nu2 <- 1/h1 + 1/h2

    # Calculate U
    Z <- diag(nh) - matrix(1/(nh), nrow=nh, ncol=nh) - 1/nu2 * crossprod(t(nu))
    Z[abs(Z) < 10^(-10)] <- 0
    U <- svd(Z)$u[, 1:(nh - 2)]
    U[abs(U) < 10^(-10)] <- 0
    
    # Construct matrix X s.t. x_t = X_t^T (1,phi,Psi)
    X <- cbind(x, matrix(0, nrow=n, ncol=nh - 1))
    X[(tau - h1 + 1):(tau + h2), 1] <- mean(x[(tau - h1 + 1):(tau + h2)]) # replace elements in (tau - h1, tau + h2) with mean
    X[(tau - h1 + 1):(tau + h2), 2] <- 1/nu2 * nu # constant of phi
    X[(tau - h1 + 1):(tau + h2), 3:ncol(X)] <- U # constants of Psi
    colnames(X) <- c("x0", "phi", paste0("psi", 1:(nh-2)))
    
    # Generate psi values
    if ( include_original ){
      if ( N > 1 ){
        Psi <- rbind(c(t(U) %*% x[(tau - h1 + 1):(tau + h2)]), matrix(rnorm((nh - 2) * (N - 1), sd=sqrt(sigma2)), nrow=(N - 1)))
      } else {
        Psi <- matrix(c(t(U) %*% x[(tau - h1 + 1):(tau + h2)]), nrow=1)
      }
    } else {
      Psi <- matrix(rnorm((nh - 2) * N, sd=sqrt(sigma2)), nrow=N)
    }
    
  } else if ( model == "var" ){

    if ( is.null(C0) ){
      C0 <- sum(x[(tau - h1 + 1):(tau + h2)]^2)
    }
    if ( is.null(C1) ){
      C1 <- sum(x[(tau - h1 + 1):tau]^2)
    }

    # Generate Psi values
    # (in the paper the W's follow a Beta distribution, but this should achieve the same result)
    if ( include_original ){
      if ( N > 1 ){
        Psi <- rbind(x[(tau - h1 + 1):(tau + h2)],
        cbind(sqrt(DirichletReg::rdirichlet(N - 1, rep(1, h1)) * C1), sqrt(DirichletReg::rdirichlet(N - 1, rep(1, h2)) * (C0 - C1))))
      } else {
        Psi <- matrix(x[(tau - h1 + 1):(tau + h2)], nrow=1)
      }
    } else {
      Psi <- cbind(sqrt(DirichletReg::rdirichlet(N, rep(1, h1)) * C1), sqrt(DirichletReg::rdirichlet(N, rep(1, h2)) * (C0 - C1)))
    }
    # (the square roots should probably be a mix of +ve and -ve, but as we square them this shouldn't matter)

  } else if ( model == "slope" ){
    
    nu <- calculate_nu_slope(n, tau - h1 + 1, tau + h2, tau, return_full=TRUE)
    nu <- nu / sqrt(sum(nu^2)) # I think sum(nu^2) should be 1, but it isn't in practice - check code?
    
    # Calculate U
    m11 <- 1/6 * ((tau + h2) * (tau + h2 + 1) * (2 * (tau + h2) + 1) - (tau - h1) * (tau - h1 + 1) * (2 * (tau - h1) + 1))
    m12 <- -1/2 * ((tau + h2) * (tau + h2 + 1) - (tau - h1) * (tau - h1 - 1))
    m22 <- h1 + h2
    M1 <- matrix(c(m11, m12, m12, m22), nrow=2) / (m22 * m11 - m12^2)
    M2 <- rbind(rep(1, nh), (tau - h1 + 1):(tau + h2))
    V <- svd(M1 %*% M2)$v
    Z <- diag(nh) - crossprod(t(V)) - crossprod(t(nu[(tau - h1 + 1):(tau + h2)]))
    Z[abs(Z) < 10^(-10)] <- 0
    U <- svd(Z)$u[, 1:(nh - 3)]
    U[abs(U) < 10^(-10)] <- 0
    
    # Construct matrix X s.t. x_t = X_t^T (1,phi,Psi)
    X <- cbind(x, matrix(0, nrow=n, ncol=nh - 2))
    X[(tau - h1 + 1):(tau + h2), 1] <- crossprod(t(V)) %*% x[(tau - h1 + 1):(tau + h2)] # constant
    X[(tau - h1 + 1):(tau + h2), 2] <- nu[(tau - h1 + 1):(tau + h2)] # coefficient of phi
    X[(tau - h1 + 1):(tau + h2), 3:(nh - 1)] <- U # coefficient of psi
    colnames(X) <- c("x0", "phi", paste0("psi", 1:(nh - 3)))
    
    # Generate psi values
    if ( include_original ){
      if ( N > 1 ){
        Psi <- rbind(c(t(U) %*% x[(tau - h1 + 1):(tau + h2)]), matrix(rnorm((nh - 3) * (N - 1), sd=sqrt(sigma2)), nrow=(N - 1)))
          # nb. check distribution is correct?
      } else {
        Psi <- matrix(c(t(U) %*% x[(tau - h1 + 1):(tau + h2)]), nrow=1)
      }
    } else {
      Psi <- matrix(rnorm((nh - 3) * N, sd=sqrt(sigma2)), nrow=N)
    }
    
  }

  if ( model %in% c("mean", "slope") ){
    return(list(X=X, U=U, Psi=Psi))
  } else if ( model == "var" ){
    return(list(X=NA, U=NA, Psi=Psi))
  }

}