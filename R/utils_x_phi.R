# Calculate X'(phi) = X - (||\nu||_2^2)^{-1} \nu \nu^T X + (||\nu||_2^2)^{-1} \nu \phi
# (for given phi value).

calculate_x_phi <- function(x, nu, phi, nu2=NULL, phi_obs=NULL){

  if ( is.null(nu2) ){
    nu2 <- sum(nu^2)
  }
  if ( is.null(phi_obs) ){
    phi_obs <- c(t(nu) %*% x)
  }
  x2 <- x + nu2^(-1) * (phi - phi_obs) * nu

  return(x2)

}


calculate_x_phi_var <- function(x, b0, phi, phi_obs, h1, h2){
  
  x2 <- x
  x2[(b0 - h1 + 1):b0] <- sqrt(phi / phi_obs) * x[(b0 - h1 + 1):b0]
  x2[(b0 + 1):(b0 + h2)] <- sqrt((1 - phi) / (1 - phi_obs)) * x[(b0 + 1):(b0 + h2)]

  return(x2)
  
}


calculate_x_phi_autocor <- function( x, nu, phi, Sigma, Sigma_nu=NULL, phi_obs=NULL, nu2=NULL ){

  if ( is.null(phi_obs) ){
    phi_obs <- c(t(nu) %*% x)
  }
  if ( is.null(nu2) ){
    nu2 <- c(t(nu) %*% Sigma %*% nu)
  }
  if ( is.null(Sigma_nu) ){
    Sigma_nu <- Sigma %*% nu
  }
  
  x2 <- x + (phi - phi_obs) / nu2 * Sigma_nu

  return(x2)

}