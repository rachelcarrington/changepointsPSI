# Binary segmentation (bs), wild binary segmentation (wbs), narrowest over threshold (not)
# change in mean, slope, variance

calculate_interval_bs <- function(x, b, d, b0=NULL, nu=NULL, model="mean", phi_var=NULL, phi_obs=NULL, threshold=0, maxiter=NULL, h1=NULL, h2=NULL, C0=NULL, C0_nu=NULL, 
                                  icss=FALSE, autocor=FALSE, Sigma=NULL ){
  # x: vector of data
  # b: vector of changepoints
  # d: vector of changepoint directions, \in (+1, -1)
  # b0: changepoint of interest, assumed to be b[1] if unspecified
  # nu: nu
  # model: "mean", "var", or "slope"
  # phi_var = sum(nu^2) (or \nu^T Sigma \nu if autocorrelation)
  # phi_obs : phi_obs
  # threshold
  # maxiter
  # h1: left window
  # h2: right window
  # C0: vector of cumulative sums of x (mean) or x^2 (var)
  # icss: TRUE for ICSS CUSUM, FALSE for standard CUSUM
  # autocor: logical, whether data is autocorrelated

  # Check that nu is supplied, if needed.
  if ( is.null(nu) ){
    if ( model == "mean" ){
      stop("nu must be specified if model == \"mean\".")
    } else if ( model == "slope" ){
      stop("nu must be specified if model == \"slope\".")
    }
  }
  
  # Autocorrelation can only be used if the model is mean (need to change this to also include slope)
#  if ( autocor & !(model == "mean") ){
#    warning("Ignoring autocorrelation.")
#    autocor <- FALSE
#  }

  n <- length(x)

  # Set maximum number of iterations
  if ( is.null(maxiter) ){
    if ( model == "mean" || model == "var" ){
      maxiter <- n - 1
    } else if ( model == "slope" ){
      maxiter <- n - 2
    }
  }
  
  # Calculate ||nu^2||
  if ( model == "mean" || model == "slope" ){
    if ( is.null(phi_var) ){
      if ( !autocor ){
        phi_var <- sum(nu^2)
      } else {
        phi_var <- c(t(nu) %*% Sigma %*% nu)
      }
    }
  }
  
  # Calculate observed value of the test statistic
  if ( model == "mean" || model == "slope" ){
    if ( is.null(phi_obs) ){
      phi_obs <- c(t(nu) %*% x)
    }
  }

  # Set changepoint of interest for variance model, if not given
  if ( model == "var" ){
    if ( is.null(b0)){
      b0 <- b[1]
    }
#    hvals <- get_h1_h2(h, b, n, b0)
#    h1 <- hvals[1]
#    h2 <- hvals[2]
  }
  
  if ( model == "slope" & autocor ){
    Sigma_nu <- Sigma %*% nu
  }


  # Calculate CUSUM statistics in terms of phi
  if ( model == "mean" ){
    cs <- cusum_phi_vec(x, nu, phi_var, phi_obs, C0=C0, C0_nu=C0_nu, icss=icss, autocor=autocor, Sigma=Sigma)
  } else if ( model == "slope" ){ 
    lambda <- calculate_nu_slope(n)
    if ( !autocor ){
      cs <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
    } else {
      cs <- t(lambda) %*% cbind(x - Sigma_nu * phi_obs / phi_var, Sigma_nu / phi_var)
    }
  } else if ( model == "var" ){
    cs <- cusum_phi_vec_var(x, b0, h1, h2, C0=C0, icss=icss)
  }

  # If no CPs (threshold must be specified):
  if ( length(b) == 0 ){

    # We need all |Ct|'s to be < threshold
    # We can ignore those for which cs[,2] = 0, since these are constant in phi
    z <- abs(cs[,2]) > 10^(-10)
    inequalities <- c((threshold - cs[z, 1]) / cs[z, 2], ((-1) * threshold - cs[z, 1]) / cs[z, 2])
    signs <- c(ifelse(cs[z, 2] > 0, -1, 1), ifelse(cs[z, 2] > 0, 1, -1))
    max_lower_bound <- max(inequalities[signs == 1])
    min_upper_bound <- min(inequalities[signs == -1])

  } else {

    # For the first CP:
    # (-d1) * C(b1) > threshold
    if ( model == "mean" || model == "var" ){
      tau <- b[1]
    } else if ( model == "slope" ){
      tau <- b[1] - 1
    }
    if ( abs(cs[tau, 2]) > 10^(-10) ){
      inequalities <- (-d[1] * threshold - cs[tau, 1]) / cs[tau, 2]
      signs <- ifelse(cs[tau, 2] > 0, -d[1], d[1])
      if ( model == "mean" || model == "slope" ){
        if ( signs == 1 ){
          max_lower_bound <- inequalities
          min_upper_bound <- Inf
        } else {
          max_lower_bound <- -Inf
          min_upper_bound <- inequalities
        }
      } else if ( model == "var" ){
        if ( signs == 1 ){
          max_lower_bound <- max(0, inequalities)
          min_upper_bound <- 1
        } else {
          max_lower_bound <- 0
          min_upper_bound <- min(1, inequalities)
        }
      }
    } else {
      if ( model == "mean" || model == "slope" ){
        max_lower_bound <- -Inf
        min_upper_bound <- Inf
      } else if ( model == "var" ){
        max_lower_bound <- 0
        min_upper_bound <- 1
      }
    }

    # |C(b1)| > +/- C(t) for t \neq b1
    inequalities <- cbind((cs[-tau, 1] - cs[tau, 1]) / (-cs[-tau, 2] + cs[tau, 2]),
                           (-1) * (cs[-tau, 1] + cs[tau, 1])/(cs[-tau, 2] + cs[tau, 2]))
    signs <- cbind(ifelse(-cs[-tau, 2] + cs[tau, 2] > 0, -d[1], d[1]),
                    ifelse(cs[-tau, 2] + cs[tau, 2] > 0, -d[1], d[1]))

    # Remove entries where we have divided by 0
    signs[abs(-cs[-tau, 2] + cs[tau, 2]) <= 10^(-10), 1] <- NA
    signs[abs(cs[-tau, 2] + cs[tau, 2]) <= 10^(-10), 2] <- NA

    # Update bounds
    max_lower_bound <- ifelse(sum(signs == 1, na.rm=TRUE) > 0, max(max_lower_bound, max(inequalities[signs == 1], na.rm=TRUE)), max_lower_bound)
    min_upper_bound <- ifelse(sum(signs == -1, na.rm=TRUE) > 0, min(min_upper_bound, min(inequalities[signs == -1], na.rm=TRUE)), min_upper_bound)

    # For later CPs, if any:
    if ( length(b) > 1 & maxiter > 1 ){

      for ( k in 2:min(length(b), maxiter) ){

        # Split data at CPs
        previous_cps <- sort(b[1:(k - 1)])
        s1 <- c(1, previous_cps + 1)
        e1 <- c(previous_cps, n)

        # Calculate CUSUM statistic in terms of phi, for each interval
        cs2 <- matrix(NA, nrow=n - 1, ncol=2)
        for ( m in 1:length(s1) ){
          if ( e1[m] > s1[m] + 1 ){  # (intervals may have length 0 if we have consecutive CPs)
            if ( model == "mean" ){
              cs2[s1[m]:(e1[m] - 1), ] <- cusum_phi_vec(x, nu, nu2=phi_var, phi_obs=phi_obs, s=s1[m], e=e1[m], C0=C0, C0_nu=C0_nu, icss=icss, autocor=autocor, Sigma=Sigma)
            } else if ( model == "slope" ){
              lambda <- calculate_nu_slope(n, s=s1[m], e=e1[m])
#              cs2[(s1[m] + 1):(e1[m] - 1), ] <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
              if ( !autocor ){
                cs2[(s1[m] + 1):(e1[m] - 1), ] <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
              } else {
                cs2[(s1[m] + 1):(e1[m] - 1), ] <- t(lambda) %*% cbind(x - Sigma_nu * phi_obs / phi_var, Sigma_nu / phi_var)
              }
            } else if ( model == "var" ){
              cs2[s1[m]:(e1[m] - 1), ] <- cusum_phi_vec_var(x, b0, h1, h2, s=s1[m], e=e1[m], C0=C0, icss=icss)
            }
          }
        }

        # We need -d[k] * C(b[k]) to be greater than 0 or the threshold
        if ( abs(cs2[b[k], 2]) > 10^(-10) ){
          if ( is.null(threshold) ){
            inequalities <- (-1) * cs2[b[k], 1] / cs2[b[k], 2]
          } else {
            inequalities <- (-d[k] * threshold - cs2[b[k], 1]) / cs2[b[k], 2]
          }
          signs <- ifelse(cs2[b[k], 2] > 0, -d[k], d[k])
          if ( signs==1 ){
            max_lower_bound <- max(max_lower_bound, inequalities)
          } else {
            min_upper_bound <- min(min_upper_bound, inequalities)
          }
        }

        # Also, need |C(b[k])| > |C(b[t])| for t \neq k
        inequalities <- cbind((cs2[-b[k], 1] - cs2[b[k], 1]) / (-cs2[-b[k], 2] + cs2[b[k], 2] ),
                               (-1)*(cs2[-b[k], 1] + cs2[b[k], 1])/(cs2[-b[k], 2] + cs2[b[k], 2]))
        signs <- cbind(ifelse(-cs2[-b[k], 2] + cs2[b[k], 2] > 0, -d[k], d[k]),
                        ifelse(cs2[-b[k], 2] + cs2[b[k], 2] > 0, -d[k], d[k]))

        # Remove entries where we have divided by 0
        signs[abs(-cs2[-b[k], 2] + cs2[b[k], 2]) <= 10^(-10), 1] <- NA
        signs[abs(cs2[-b[k], 2] + cs2[b[k], 2]) <= 10^(-10), 2] <- NA

        # Remove NAs & recalculate bounds
        max_lower_bound <- ifelse(sum(signs == 1, na.rm=TRUE) > 0, max(max_lower_bound, max(inequalities[signs == 1], na.rm=TRUE)), max_lower_bound)
        min_upper_bound <- ifelse(sum(signs == -1, na.rm=TRUE) > 0, min(min_upper_bound, min(inequalities[signs == -1], na.rm=TRUE)), min_upper_bound)

      }
    }
    
    # Check that the next CP is below the threshold, if applicable
    if ( maxiter > length(b) ){
      
      # Check that for the next k, all |Ct|'s are below the threshold
      
      # Split data at CPs
      previous_cps <- sort(b)
      s1 <- c(1, previous_cps + 1)
      e1 <- c(previous_cps, n)
      
      # Calculate CUSUM statistic in terms of phi, for each interval
      cs2 <- matrix(NA, nrow=n - 1, ncol=2)
      for ( m in 1:length(s1) ){
        if ( e1[m] > s1[m] + 1 ){  # (intervals may have length 0 if we have consecutive CPs)
          if ( model == "mean" ){
            cs2[s1[m]:(e1[m] - 1), ] <- cusum_phi_vec(x, nu, nu2=phi_var, phi_obs=phi_obs, s=s1[m], e=e1[m], C0=C0, C0_nu=C0_nu, icss=icss, autocor=autocor, Sigma=Sigma)
          } else if ( model == "slope" ){
            lambda <- calculate_nu_slope(n, s=s1[m], e=e1[m])
            if ( !autocor ){
              cs2[(s1[m] + 1):(e1[m] - 1), ] <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
            } else {
              cs2[(s1[m] + 1):(e1[m] - 1), ] <- t(lambda) %*% cbind(x - Sigma_nu * phi_obs / phi_var, Sigma_nu / phi_var)
            }
#            cs2[(s1[m] + 1):(e1[m] - 1), ] <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
          } else if ( model == "var" ){
            cs2[s1[m]:(e1[m] - 1), ] <- cusum_phi_vec_var(x, b0, h1, h2, s=s1[m], e=e1[m], C0=C0, icss=icss)
          }
        }
      }
      
      z <- abs(cs2[,2]) > 10^(-10)
      inequalities <- c((threshold - cs2[z, 1]) / cs2[z, 2], ((-1) * threshold - cs2[z, 1]) / cs2[z, 2])
      signs <- c(ifelse(cs2[z, 2] > 0, -1, 1), ifelse(cs2[z, 2] > 0, 1, -1))
      max_lower_bound <- ifelse(sum(signs == 1, na.rm=TRUE) > 0, max(c(max_lower_bound, inequalities[signs == 1]), na.rm=TRUE), max_lower_bound)
      min_upper_bound <- ifelse(sum(signs == -1, na.rm=TRUE) > 0, min(c(min_upper_bound, inequalities[signs == -1]), na.rm=TRUE), min_upper_bound)
      
    }

  }

  return(c(max_lower_bound, min_upper_bound))
}




calculate_interval_wbs <- function(x, b, d, s, e, random_samples, nu=NULL, phi_var=NULL, phi_obs=NULL, threshold=0, maxiter=NULL, tau0=NULL, model="mean", h1=NULL, h2=NULL,
                                   C0=NULL, C0_nu=NULL, icss=FALSE, autocor=FALSE, Sigma=NULL ){
  
  n <- length(x)

  if ( is.null(nu) ){
    if ( model == "mean" ){
      stop("nu must be specified if model == \"mean\".")
    } else if ( model == "slope" ){
      stop("nu must be specified if model == \"slope\".")
    }
  }

  if ( model == "var" & is.null(tau0) ){
    tau0 <- b[1]
  }
  
  if ( model == "mean" || model == "slope" ){
    if ( is.null(phi_var) ){
      if ( !autocor ){
        phi_var <- sum(nu^2)
      } else {
        phi_var <- c(t(nu) %*% Sigma %*% nu)
      }
    }
  }
  
  if ( model == "mean" || model == "slope" ){
    if ( is.null(phi_obs) ){
      phi_obs <- c(t(nu) %*% x)
    }
  }

  max_lower_bound <- ifelse(model == "var", 0, -Inf)
  min_upper_bound <- ifelse(model == "var", 1, Inf)
  
  # If no CPs:
  if ( sum(!is.na(b)) == 0 ){

    # For all intervals (s,e), we need |C_{s,e}(t)| < \lambda for all t
    for ( ind in 1:nrow(random_samples) ){
      if ( model == "mean" ){
        cs <- cusum_phi_vec(x, nu, phi_var, phi_obs, s=random_samples[ind, 1], e=random_samples[ind, 2], C0=C0, C0_nu=C0_nu, icss=icss, autocor=autocor, Sigma=Sigma)
      } else if ( model == "slope" ){
        lambda <- calculate_nu_slope(n, s=random_samples[ind, 1], e=random_samples[ind, 2])
#        cs <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
        if ( !autocor ){
          cs <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
        } else {
          cs <- t(lambda) %*% cbind(x - Sigma_nu * phi_obs / phi_var, Sigma_nu / phi_var)
        }
      } else if ( model == "var" ){
        cs <- cusum_phi_vec_var(x, tau0, h1, h2, s=random_samples[ind, 1], e=random_samples[ind, 2], C0=C0, icss=icss)
      }
      z <- abs(cs[,2]) > 10^(-10) # if cs[i,2] = 0, then this inequality is constant in phi
      inequalities <- c((threshold - cs[z, 1]) / cs[z, 2], (-threshold - cs[z,1]) / cs[z,2])
      signs <- c(ifelse(cs[z, 2] > 0, -1, 1), ifelse(cs[z, 2] > 0, 1, -1))
      max_lower_bound <- ifelse(sum(signs == 1, na.rm=TRUE) > 0, max(c(max_lower_bound, inequalities[signs == 1]), na.rm=TRUE), max_lower_bound)
      min_upper_bound <- ifelse(sum(signs == -1, na.rm=TRUE) > 0, min(c(min_upper_bound, inequalities[signs == -1]), na.rm=TRUE), min_upper_bound)
    }

  } else {

    for ( k in 1:min(length(b), maxiter) ){
       
      # Remove intervals containing previous changepoint
      if ( k > 1 ){
        random_samples <- random_samples[!(random_samples[,1] <= b[k - 1] & random_samples[,2] > b[k - 1]), ]
      }
      
      # Also remove CP-containing interval (we can consider this separately from the others)
      random_samples <- random_samples[!(random_samples[,1] == s[k] & random_samples[,2] == e[k]), ]
    
      # For interval containing CP, we need |C(t)| < |C(b1)| for t \neq b1
      if ( model == "mean" ){
        cs <- cusum_phi_vec(x, nu, phi_var, phi_obs, s=s[k], e=e[k], C0=C0, C0_nu=C0_nu, icss=icss, autocor=autocor, Sigma=Sigma)
      } else if ( model == "slope" ){
#        lambda <- create_lambda_change_in_slope(n, s=s[k], e=e[k])
        lambda <- calculate_nu_slope(n, s=s[k], e=e[k])
        cs <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
      } else if ( model == "var" ){
        cs <- cusum_phi_vec_var(x, tau0, h1, h2, s=s[k], e=e[k], C0=C0, icss=icss)
      }

      if ( model == "mean" || model == "var" ){
        tau <- b[k] - s[k] + 1
      } else if ( model == "slope" ){
        tau <- b[k] - s[k]
      }
    
      # (-d[k]) * C(b[k]) > threshold
      if ( abs(cs[tau, 2]) > 10^(-10) ){
        inequalities <- (-d[k] * threshold - cs[tau, 1]) / cs[tau, 2]
        signs <- ifelse(cs[tau, 2] > 0, -d[k], d[k])
        if ( signs == 1 ){
          max_lower_bound <- max(max_lower_bound, inequalities)
        } else {
          min_upper_bound <- min(min_upper_bound, inequalities)
        }
      }
    
      # | C(b1) | > +/- C(t) for t \neq b[k], on interval (s[k], e[k])
      inequalities <- cbind((cs[-tau, 1] - cs[tau, 1])/(-cs[-tau, 2] + cs[tau, 2]), 
                             (-1) * (cs[-tau, 1] + cs[tau, 1])/(cs[-tau, 2] + cs[tau, 2]))
      signs <- cbind(ifelse(-cs[-tau, 2] + cs[tau, 2] > 0, -d[k], d[k]),
                      ifelse(cs[-tau, 2] + cs[tau, 2] > 0, -d[k], d[k]))
      
      # Remove entries where we have divided by 0
      signs[abs(-cs[-tau, 2] + cs[tau, 2]) <= 10^(-10), 1] <- NA
      signs[abs(cs[-tau, 2] + cs[tau, 2]) <= 10^(-10), 2] <- NA
    
      # Update bounds
      max_lower_bound <- ifelse(sum(signs == 1, na.rm=TRUE) > 0, max(c(max_lower_bound, inequalities[signs == 1]), na.rm=TRUE), max_lower_bound)
      min_upper_bound <- ifelse(sum(signs == -1, na.rm=TRUE) > 0, min(c(min_upper_bound, inequalities[signs == -1]), na.rm=TRUE), min_upper_bound)
    
      # Check other intervals: | C_{s1, e1}(b1) | > |C_{s,e} (t)| for all s, e, t
      for ( ind in 1:nrow(random_samples) ){
        s2 <- random_samples[ind, 1]
        e2 <- random_samples[ind, 2]
        if ( model == "mean" ){
          cs2 <- cusum_phi_vec(x, nu, phi_var, phi_obs, s=s2, e=e2, C0=C0, C0_nu=C0_nu, icss=icss, autocor=autocor, Sigma=Sigma)
        } else if ( model == "slope" ){
          lambda <- calculate_nu_slope(n, s=s2, e=e2)
#          cs2 <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
          if ( !autocor ){
            cs2 <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
          } else {
            cs2 <- t(lambda) %*% cbind(x - Sigma_nu * phi_obs / phi_var, Sigma_nu / phi_var)
          }
        } else if ( model == "var" ){
          cs2 <- cusum_phi_vec_var(x, tau0, h1, h2, s=s2, e=e2, C0=C0, icss=icss)
        }
 
        inequalities <- cbind((cs2[,1] - cs[tau,1])/(-cs2[,2] + cs[tau,2] ), 
                               (-1)*(cs2[,1] + cs[tau,1])/(cs2[,2] + cs[tau,2]))
        signs <- cbind(ifelse(-cs2[,2] + cs[tau,2] > 0, -d[k], d[k]),
                        ifelse(cs2[,2] + cs[tau,2] > 0, -d[k], d[k]))
        
        # Remove entries where we have divided by 0
        signs[abs(-cs2[,2] + cs[tau,2]) <= 10^(-10), 1] <- NA
        signs[abs(cs2[,2] + cs[tau,2]) <= 10^(-10), 2] <- NA
        
        # Update bounds
        max_lower_bound <- ifelse(sum(signs == 1, na.rm=TRUE) > 0, max(c(max_lower_bound, inequalities[signs == 1]), na.rm=TRUE), max_lower_bound)
        min_upper_bound <- ifelse(sum(signs == -1, na.rm=TRUE) > 0, min(c(min_upper_bound, inequalities[signs == -1]), na.rm=TRUE), min_upper_bound)
      }
 
    }
        
    if ( maxiter > length(b) ){

      # Check that for the next k, all |C_t|'s are below the threshold
      k <- length(b) + 1

      # Remove intervals which contain the previous changepoint
      random_samples <- random_samples[!(random_samples[,1] <= b[k - 1] & random_samples[,2] > b[k - 1]), , drop=FALSE]
 
      # Calculate bounds for which all C(t)'s in each RI will be below the threshold
      if ( nrow(random_samples) >= 1 ){
        for ( ind in 1:nrow(random_samples) ){
          if ( model == "mean" ){
            cs <- cusum_phi_vec(x, nu, phi_var, phi_obs, s=random_samples[ind, 1], e=random_samples[ind, 2], C0=C0, C0_nu=C0_nu, icss=icss, autocor=autocor, Sigma=Sigma)
          } else if ( model == "slope" ){
            lambda <- calculate_nu_slope(n, s=random_samples[ind, 1], e=random_samples[ind, 2])
            if ( !autocor ){
              cs <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
            } else {
              cs <- t(lambda) %*% cbind(x - Sigma_nu * phi_obs / phi_var, Sigma_nu / phi_var)
            }
#            cs <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
          } else if ( model == "var" ){
            cs <- cusum_phi_vec_var(x, tau0, h1, h2, s=random_samples[ind, 1], e=random_samples[ind, 2], C0=C0, icss=icss)
          }
          z <- abs(cs[,2]) > 10^(-10) # if cs[i,2] = 0, then this inequality is constant in phi
          inequalities <- c((threshold - cs[z,1]) / cs[z,2], (-threshold - cs[z, 1]) / cs[z, 2])
          signs <- c(ifelse(cs[z,2] > 0, -1, 1), ifelse(cs[z,2] > 0, 1, -1))
          max_lower_bound <- ifelse(sum(signs == 1, na.rm=TRUE) > 0, max(c(max_lower_bound, inequalities[signs == 1]), na.rm=TRUE), max_lower_bound)
          min_upper_bound <- ifelse(sum(signs == -1, na.rm=TRUE) > 0, min(c(min_upper_bound, inequalities[signs == -1]), na.rm=TRUE), min_upper_bound)
        }
      }

    }

  }
    
  # Inequalities are all satisfied by phi s.t. max_lower_bound < phi < min_upper_bound
  return( c(max_lower_bound, min_upper_bound) )
}





calculate_interval_not <- function(x, b, d, s, e, random_samples, threshold, maxiter, nu=NULL, phi_var=NULL, phi_obs=NULL, h1=NULL, h2=NULL, tau0=NULL, model="mean", C0=NULL,
                                   C0_nu=NULL, icss=FALSE, autocor=FALSE, Sigma=NULL){
  
  # Find values of phi which satisfy the required inequalities so that the NOT (change in mean) algorithm returns given results
  
  if ( is.null(nu) ){
    if ( model == "mean" ){
      stop("nu must be specified if model == \"mean\".")
    } else if ( model == "slope" ){
      stop("nu must be specified if model == \"slope\".")
    }
  }

  n <- length(x)

  if ( model == "mean" || model == "slope" ){
    if ( is.null(phi_var) ){
      if ( !autocor ){
        phi_var <- sum(nu^2)
      } else {
        phi_var <- c(t(nu) %*% Sigma %*% nu)
      }
    }
  }
  
  if ( model == "mean" || model == "slope" ){
    if ( is.null(phi_obs) ){
      phi_obs <- c(t(nu) %*% x)
    }
  }

  if ( model == "var" & is.null(tau0) ){
    tau0 <- b[1]
  }

  max_lower_bound <- ifelse(model == "var", 0, -Inf)
  min_upper_bound <- ifelse(model == "var", 1, Inf)
  
  if ( length(b) == 0 ){
    
    # If there are no CPs detected, we require |C_(s,e) (t)| < lambda for all (s,e) and all t
    inequalities <- signs <- numeric(0)
    for ( j in 1:nrow(random_samples) ){
      if ( model == "mean" ){
        cs <- cusum_phi_vec(x, nu, nu2=phi_var, phi_obs=phi_obs, s=random_samples[j,1], e=random_samples[j,2], C0=C0, C0_nu=C0_nu, icss=icss, autocor=autocor, Sigma=Sigma)
      } else if ( model == "slope" ){
#        lambda <- create_lambda_change_in_slope(n, s=random_samples[j,1], e=random_samples[j,2])
        lambda <- calculate_nu_slope(n, s=random_samples[j,1], e=random_samples[j,2])
#        cs <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
        if ( !autocor ){
          cs <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
        } else {
          cs <- t(lambda) %*% cbind(x - Sigma_nu * phi_obs / phi_var, Sigma_nu / phi_var)
        }
      } else if ( model == "var" ){
        cs <- cusum_phi_vec_var(x, tau0, h1, h2, s=random_samples[ind,1], e=random_samples[ind,2], C0=C0, icss=icss)
      }
      cs <- cs[abs(cs[,2]) > 10^(-10),,drop=FALSE]
      inequalities <- c((threshold - cs[,1])/cs[,2], ((-1) * threshold - cs[,1])/cs[,2])
      signs <- c(ifelse(cs[,2] > 0, -1, 1), ifelse(cs[,2] > 0, 1, -1))
      max_lower_bound <- ifelse(sum(signs == 1, na.rm=TRUE) > 0, max(c(max_lower_bound, inequalities[signs == 1]), na.rm=TRUE), max_lower_bound)
      min_upper_bound <- ifelse(sum(signs == -1, na.rm=TRUE) > 0, min(c(min_upper_bound, inequalities[signs == -1]), na.rm=TRUE), min_upper_bound)
    }
    
  } else {
    
    for ( k in 1:min(length(b), maxiter) ){
      
      # Remove intervals containing previous changepoint
      # also intervals which are smaller than previous CP-containing changepoint, as we have already checked that these are below the threshold
      if ( k > 1 ){
        random_samples <- random_samples[!(random_samples[,1] <= b[k - 1] & random_samples[,2] > b[k - 1]), ]
        random_samples <- random_samples[random_samples[,2] - random_samples[,1] >= e[k - 1] - s[k - 1], ]
      }
      
      # For narrower intervals, need |C(t)| < threshold for all t in interval
      widths <- random_samples[,2] - random_samples[,1] + 1
      cp_width <- e[k] - s[k] + 1
      smaller_intervals <-random_samples[widths < cp_width,,drop=FALSE]
      
      inequalities <- signs <- numeric(0)
      if ( nrow(smaller_intervals) >= 1 ){
        for ( j in 1:nrow(smaller_intervals) ){
          s1 <- smaller_intervals[j,1]
          e1 <- smaller_intervals[j,2]
          
          if ( model == "mean" ){
            cs <- cusum_phi_vec(x, nu, phi_var, phi_obs, s=s1, e=e1, C0=C0, C0_nu=C0_nu, icss=icss, autocor=autocor, Sigma=Sigma)
          } else if ( model == "slope" ){
#            lambda <- create_lambda_change_in_slope(n, s=s1, e=e1)
            lambda <- calculate_nu_slope(n, s=s1, e=e1)
#            cs <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
            if ( !autocor ){
              cs <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
            } else {
              cs <- t(lambda) %*% cbind(x - Sigma_nu * phi_obs / phi_var, Sigma_nu / phi_var)
            }
          }  else if ( model == "var" ){
            cs <- cusum_phi_vec_var(x, tau0, h1, h2, s=s1, e=e1, C0=C0, icss=icss)
          }
          cs <- cs[abs(cs[,2]) > 10^(-10),,drop=FALSE] # if cs[t,2] = 0 then C(t) is constant in phi
          inequalities <- c(inequalities, (threshold - cs[,1]) / cs[,2], ((-1) * threshold - cs[,1]) / cs[,2])
          signs <- c(signs, ifelse(cs[,2] > 0, -1, 1), ifelse(cs[,2] > 0, 1, -1))
        }
        max_lower_bound <- ifelse(sum(signs == 1, na.rm=TRUE) > 0, max(c(max_lower_bound, inequalities[signs == 1]), na.rm=TRUE), max_lower_bound)
        min_upper_bound <- ifelse(sum(signs == -1, na.rm=TRUE) > 0, min(c(min_upper_bound, inequalities[signs == -1]), na.rm=TRUE), min_upper_bound)
      }
      
      # For interval containing changepoint, we need |C(t)| < C(b1) for t \neq b1
      if ( model == "mean" ){
        cs <- cusum_phi_vec(x, nu, nu2=phi_var, phi_obs=phi_obs, s=s[k], e=e[k], C0=C0, C0_nu=C0_nu, icss=icss, autocor=autocor, Sigma=Sigma)
      } else if ( model == "slope" ){
#        lambda <- create_lambda_change_in_slope(n, s=s[k], e=e[k])
        lambda <- calculate_nu_slope(n, s=s[k], e=e[k])
#        cs <- rbind(c(NA, NA), t(lambda) %*% cbind(x - nu * phi_obs, nu))
        if ( !autocor ){
          cs <- rbind(c(NA, NA), t(lambda) %*% cbind(x - nu * phi_obs, nu))
        } else {
          cs <- rbind(c(NA, NA), t(lambda) %*% cbind(x - Sigma_nu * phi_obs / phi_var, Sigma_nu / phi_var))
        }
      }  else if ( model == "var" ){
        cs <- cusum_phi_vec_var(x, tau0, h1, h2, s=s[k], e=e[k], C0=C0, icss=icss)
      }
      
      # C(b[k]) > (-d1) * threshold
      if ( abs(cs[b[k] - s[k] + 1, 2]) > 10^(-10) ){ # if denominator is 0, inequality is always satisfied so we can ignore it
        if ( d[k] == 1 ){
          inequalities <- (-1) * (threshold + cs[b[k] - s[k] + 1, 1])/cs[b[k] - s[k] + 1, 2]
        } else {
          inequalities <- (threshold - cs[b[k] - s[k] + 1, 1])/cs[b[k] - s[k] + 1, 2]
        }
        signs <- ifelse(cs[b[k] - s[k] + 1, 2] > 0, -d[k], d[k])
        max_lower_bound <- ifelse(sum(signs == 1, na.rm=TRUE) > 0, max(c(max_lower_bound, inequalities[signs == 1]), na.rm=TRUE), max_lower_bound)
        min_upper_bound <- ifelse(sum(signs == -1, na.rm=TRUE) > 0, min(c(min_upper_bound, inequalities[signs == -1]), na.rm=TRUE), min_upper_bound)
      }
      
      # |C(b[k])| > +/- C(t) for t \neq b[k]
      inequalities <- cbind((cs[,1] - cs[b[k] - s[k] + 1, 1])/(cs[b[k] - s[k] + 1, 2] - cs[,2]), 
                            (-1)*(cs[,1] + cs[b[k] - s[k] + 1, 1])/(cs[b[k] - s[k] + 1, 2] + cs[,2]))
      signs <- cbind(ifelse(cs[b[k] - s[k] + 1, 2] - cs[,2] > 0, -d[k], d[k]), ifelse(cs[b[k] - s[k] + 1, 2] + cs[,2] > 0, -d[k], d[k]))
      
      # Remove b[k] since we don't need this (will probably return NaN or Inf otherwise)
      signs[b[k] - s[k] + 1, ] <- NA
      signs[abs(cs[,2] - cs[b[k] - s[k] + 1, 2]) <= 10^(-10)] <- NA
      signs[abs(cs[,2] + cs[b[k] - s[k] + 1, 2]) <= 10^(-10)] <- NA
      max_lower_bound <- ifelse(sum(signs == 1, na.rm=TRUE) > 0, max(c(max_lower_bound, inequalities[signs == 1]), na.rm=TRUE), max_lower_bound)
      min_upper_bound <- ifelse(sum(signs == -1, na.rm=TRUE) > 0, min(c(min_upper_bound, inequalities[signs == -1]), na.rm=TRUE), min_upper_bound)
      
      # For intervals of the same width, we need |C(t)| < C_(s,e) (b1) for all t in interval
      same_intervals <- random_samples[widths == cp_width,,drop=FALSE]
      same_intervals <- same_intervals[!(same_intervals[,1] == s[k] & same_intervals[,2] == e[k]),,drop=FALSE]
      if ( nrow(same_intervals) >= 1 ){
        for ( j in nrow(same_intervals) ){
          if ( model == "mean" ){
            cs2 <- cusum_phi_vec(x, nu, phi_var, phi_obs, s=same_intervals[j, 1], e=same_intervals[j, 2], C0=C0, C0_nu=C0_nu, icss=icss, autocor=autocor, Sigma=Sigma)
          } else if ( model == "slope" ){
#            lambda <- create_lambda_change_in_slope(n, s=same_intervals[j, 1], e=same_intervals[j, 2])
            lambda <- calculate_nu_slope(n, s=same_intervals[j, 1], e=same_intervals[j, 2])
#            cs2 <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
            if ( !autocor ){
              cs2 <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
            } else {
              cs2 <- t(lambda) %*% cbind(x - Sigma_nu * phi_obs / phi_var, Sigma_nu / phi_var)
            }
          } else if ( model == "var" ){
            cs2 <- cusum_phi_vec_var(x, tau0, h1, h2, s=same_intervals[j, 1], e=same_intervals[j, 2], C0=C0, icss=icss)
          }
          cs2 <- cs2[abs(cs2[,2]) >= 10^(-10),,drop=FALSE]
          if ( nrow(cs2) > 0 ){
            inequalities <- cbind((cs2[,1] - cs[b[k] - s[k] + 1, 1]) / (cs[b[k] - s[k] + 1, 2] - cs2[,2]), 
                              (-1)*(cs2[,1] + cs[b[k] - s[k] + 1, 1])/(cs[b[k] - s[k] + 1, 2] + cs2[,2]))
            signs <- cbind(ifelse(cs[b[k] - s[k] + 1, 2] - cs2[,2] > 0, -d[k], d[k]), ifelse(cs[b[k] - s[k] + 1, 2] + cs2[,2] > 0, -d[k], d[k]))
            
            # Remove entries where we have divided by 0
            signs[abs(cs2[,2] - cs[b[k] - s[k] + 1, 2]) <= 10^(-10), 1] <- NA
            signs[abs(cs2[,2] + cs[b[k] - s[k] + 1, 2]) <= 10^(-10), 2] <- NA
            
            max_lower_bound <- ifelse(sum(signs == 1, na.rm=TRUE) > 0, max(c(max_lower_bound, inequalities[signs == 1]), na.rm=TRUE), max_lower_bound)
            min_upper_bound <- ifelse(sum(signs == -1, na.rm=TRUE) > 0, min(c(min_upper_bound, inequalities[signs == -1]), na.rm=TRUE), min_upper_bound)
          }
        }
      }
      
    }
    
    # Make sure no more changepoints are detected
    if ( maxiter > length(b) ){
      
      # Remove intervals containing last changepoint
      k <- length(b) + 1
      random_samples <- random_samples[!(random_samples[,1] <= b[k - 1] & random_samples[,2] > b[k - 1]),, drop=FALSE]
      
      # Check that all remaining intervals are below the threshold
      if ( nrow(random_samples) >= 1 ){
        inequalities <- signs <- numeric(0)
        for ( j in 1:nrow(random_samples) ){  
          if ( model == "mean" ){
            cs <- cusum_phi_vec(x, nu, nu2=phi_var, phi_obs=phi_obs, s=random_samples[j,1], e=random_samples[j,2], C0=C0, C0_nu=C0_nu, icss=icss, autocor=autocor, Sigma=Sigma)
          } else if ( model == "slope" ){
#            lambda <- create_lambda_change_in_slope(n, s=random_samples[j, 1], e=random_samples[j, 2])
            lambda <- calculate_nu_slope(n, s=random_samples[j, 1], e=random_samples[j, 2])
#            cs <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
            if ( !autocor ){
              cs <- t(lambda) %*% cbind(x - nu * phi_obs, nu)
            } else {
              cs <- t(lambda) %*% cbind(x - Sigma_nu * phi_obs / phi_var, Sigma_nu / phi_var)
            }
          } else if ( model == "var" ){
            cs <- cusum_phi_vec_var(x, tau0, h1, h2, s=random_samples[j, 1], e=random_samples[j, 2], C0=C0, icss=icss)
          }
          cs <- cs[abs(cs[,2]) > 10^(-10),,drop=FALSE]
          inequalities <- c(inequalities, (threshold - cs[,1]) / cs[,2], ((-1) * threshold - cs[,1]) / cs[,2])
          signs <- c(signs, ifelse(cs[,2] > 0, -1, 1), ifelse(cs[,2] > 0, 1, -1))
        }
        max_lower_bound <- ifelse(sum(signs == 1, na.rm=TRUE) > 0, max(c(max_lower_bound, inequalities[signs == 1]), na.rm=TRUE), max_lower_bound)
        min_upper_bound <- ifelse(sum(signs == -1, na.rm=TRUE) > 0, min(c(min_upper_bound, inequalities[signs == -1]), na.rm=TRUE), min_upper_bound)
      }
      
    }
    
  }
  
  
  return(c(max_lower_bound, min_upper_bound))
}