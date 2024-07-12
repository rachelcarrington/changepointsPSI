#' Calculate S
#'
#' @description For a given dataset and changepoint algorithm (one of binary segmentation, wild binary segmentation, or narrowest over threshold),
#' calculate the possible outcomes that can be obtained by applying the algorithm to
#' the perturbed dataset \eqn{X'(\phi)}, and the intervals of \eqn{\phi} values that lead to each possible outcome. See below for more details.
#'
#' @param x Vector of data.
#' @param results Output of \code{find_changepoints}.
#' @param tau Changepoint of interest. Defaults to the first detected changepoint, if not given.
#' @param nu Vector that defines the test statistic, for \code{model = "mean"} or \code{model = "slope"}. See below.
#' @param nu2 Value of \eqn{||\nu||_2^2}; optional.
#' @param h Window size; only needed if \code{model = "var"}. If the window is asymmetric, this should be a vector of length 2 containing the left and right window sizes.
#' Otherwise, it should be a single integer value.
#' @param phi_obs Observed value of the test statistic; optional.
#' @param first_cp_only Logical. If \code{TRUE}, condition on the fact that the changepoint of interest is in the model; 
#' if \code{FALSE}, condition on all changepoints. Defaults to \code{FALSE}.
#'
#' @return A dataframe containing intervals with the changepoints obtained when \eqn{\phi} is in each interval.
#'
#' @details There are three options for the changepoint model:
#' \itemize{
#' \item Change in mean (\code{model = "mean"}): \eqn{X_t = \mu_t + \epsilon_t} where \eqn{\mu_t} is piecewise constant and \eqn{\epsilon_t \sim N(0, \sigma^2)}.
#' \item Change in slope (\code{model = "slope"}): same as above but with \eqn{\mu_t} piecewise linear.
#' \item Change in variance (\code{model = "var"}): \eqn{\mu_t} is assumed to be 0 for all \eqn{t}, and \eqn{\epsilon_t \sim N(0, \sigma_t^2)} where \eqn{\sigma_t}
#' is piecewise constant.
#' }
#'
#' For the change in mean (\code{model = "mean"}) and change in slope (\code{model = "slope"}), \eqn{X'(\phi)} is defined as
#' \deqn{X'(\phi) = X_{obs} + \frac{1}{||\nu||^2} \nu (\phi - \phi_{obs}).}
#' The vector \code{nu} must be supplied to the function.
#' For the change in variance model, \eqn{X'(\phi)} is defined by
#' \deqn{
#' X'(\phi) = \begin{cases}
#' X_t & t \leq \hat{\tau} - h_1 \text{ or } t \geq \hat{\tau} + h_2 + 1 \\
#' \sqrt{\frac{\phi}{\phi_{obs}}} X_t & \hat{\tau} - h_1 + 1 \leq t \leq \hat{\tau} \\
#' \sqrt{\frac{1 - \phi}{1 - \phi_{obs}}} X_t & \hat{\tau} + 1 \leq t \leq \hat{\tau} + h_2.
#' \end{cases}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(100)
#' x <- rnorm(100) + c(rep(1,50), rep(-1,50))
#' results <- find_changepoints(x, "bs")
#' b <- results$changepoints
#' print(b)
#' h <- 10
#' nu <- c(rep(0, b[1] - h), rep(1/h, h), rep(-1/h, h), rep(0, length(x) - b[1] - h))
#' calculate_S(x, results, nu=nu, first_cp_only=TRUE)
#'
#'
calculate_S <- function(x, results, tau=NULL, nu=NULL, nu2=NULL, h=NULL, phi_obs=NULL, first_cp_only=FALSE){

  # h is only needed if model == "var" -- then nu is not provided

  model <- results$params$model
  method <- results$params$method
  cp_params <- within(results$params, rm(model, method))

  if ( is.null(nu) ){
    if ( model == "mean" ){
      stop("nu must be specified if model == \"mean\".")
    } else if ( model == "slope" ){
      stop("nu must be specified if model == \"slope\".")
    }
  }

  eps0 <- ifelse(model == "var", 0.001, 0.01)

  n <- length(x)

  if ( model =="mean" & is.null(nu2) ){
    nu2 <- sum(nu^2)
  }

  if ( model %in% c("mean", "slope") & is.null(phi_obs) ){
    phi_obs <- c(t(nu) %*% x)
  }

  b <- results$changepoints
  d <- results$results$d
  maxiter <- results$params$maxiter

  if ( is.null(tau) ){ # tau is the CP of interest
    tau <- b[1]
  }

  if ( model == "var" ){
    if ( length(h) == 1 ){
      h1 <- h2 <- h
    } else if ( length(h) == 2 ){
      h1 <- h[1]
      h2 <- h[2]
    } else {
      stop("Invalid h.")
    }
    if ( is.null(phi_obs) ){
      C0 <- sum(x[(tau - h1 + 1):(tau + h2)]^2)
      C1 <- sum(x[(tau - h1 + 1):tau]^2)
      phi_obs <- C1 / C0
    }
  }


  if ( first_cp_only ){
    
    # Calculate first interval
    interval <- calculate_interval(results, nu=nu, nu2=nu2, phi_obs=phi_obs, h1=h1, h2=h2, maxiter=1, b0=tau)
    S <- matrix(c(interval, as.matrix(c(b[1], d[1]))), nrow=1)
    colnames(S) <- c("lower_lim", "upper_lim", "b1", "d1")

    # Calculate higher intervals
    max_cps <- 1
    eps <- eps0
    upper_limit <- ifelse(model == "var", 1, Inf)
    while ( max(S[,"upper_lim"]) < upper_limit ){
      if ( model == "mean" ){
        phi <- max(S[,"upper_lim"]) + eps
        x_phi <- calculate_x_phi(x, nu, phi, nu2=nu2, phi_obs=phi_obs)
      } else if ( model == "slope" ){
        phi <- max(S[,"upper_lim"]) + eps
        x_phi <- calculate_x_phi(x, nu, phi, phi_obs=phi_obs)
      } else if ( model == "var" ){
        phi <- min(max(S[,"upper_lim"]) + eps, 1)
        x_phi <- calculate_x_phi_var(x, tau, phi, phi_obs, h1, h2)
      }

      results_phi <- find_changepoints(x_phi, method, cp_params, model)
      b2 <- results_phi$changepoints
            
      if ( b[1] %in% b2 ){ # then we only have to go as far as b[1] in calculating CPs
        n.cp <- (1:length(b2))[b2 == b[1]]
      } else {
        n.cp <- maxiter
      }

      interval <- calculate_interval(results_phi, x=x, nu=nu, nu2=nu2, phi_obs=phi_obs, h1=h1, h2=h2, maxiter=n.cp, b0=tau)

      # Check this interval is the next one
      ncps_found <- min(n.cp, sum(!is.na(results_phi$changepoints)))
      if ( abs(interval[1] - max(S[,"upper_lim"])) < 10^(-10) ){
        if ( ncps_found == 0 ){
          S <- rbind(S, c(interval, rep(NA, ncol(S) - 2)))
        } else if ( ncps_found == max_cps ){
          S <- rbind(S, c(interval, b2[1:ncps_found], results_phi$results$d[1:ncps_found]))
        } else if ( ncps_found <= max_cps ){
          S <- rbind(S, c(interval, b2[1:ncps_found], rep(NA, max_cps - ncps_found), results_phi$results$d[1:ncps_found], rep(NA, max_cps - ncps_found)) )
        } else {
          S <- cbind(S[, 1:(2 + max_cps), drop=FALSE], matrix(NA, nrow=nrow(S), ncol=ncps_found - max_cps), S[, -(1:(2 + max_cps)), drop=FALSE],
                 matrix(NA, nrow=nrow(S), ncol=ncps_found - max_cps))
          colnames(S) <- c("lower_lim", "upper_lim", paste0("b", 1:((ncol(S) - 2)/2)), paste0("d", 1:((ncol(S) - 2)/2)))
          S <- rbind(S, c(interval, b2[1:ncps_found], results_phi$results$d[1:ncps_found]))
          max_cps <- ncps_found
        }
        eps <- eps0
      } else { # otherwise, try again with a smaller epsilon
        eps <- eps/2
      }

    }

    eps <- eps0
    lower_limit <- ifelse(model == "var", 0, -Inf)
    while ( min(S[,"lower_lim"]) > lower_limit ){
      if ( model == "mean" ){
        phi <- min(S[,"lower_lim"]) - eps
        x_phi <- calculate_x_phi(x, nu, phi, nu2=nu2, phi_obs=phi_obs)
      } else if ( model == "slope" ){
        phi <- max(S[,"upper_lim"]) + eps
        x_phi <- calculate_x_phi(x, nu, phi, phi_obs=phi_obs)
      } else if ( model == "var" ){
        phi <- max(min(S[,"lower_lim"]) - eps, 0)
        x_phi <- calculate_x_phi_var(x, tau, phi, phi_obs, h1, h2)
      }
      results_phi <- find_changepoints(x_phi, method, cp_params, model)
      b2 <- results_phi$changepoints      
      
      if ( b[1] %in% b2 ){ # then we only have to go as far as b[1] in calculating CPs
        n.cp <- (1:length(b2))[b2 == b[1]]
      } else {
        n.cp <- maxiter
      }

      interval <- calculate_interval(results_phi, x=x, nu=nu, nu2=nu2, phi_obs=phi_obs, h1=h1, h2=h2, maxiter=n.cp, b0=tau)

      # Check this interval is the next one
      ncps_found <- min(n.cp, sum(!is.na(b2)))
      if ( abs(interval[2] - min(S[,"lower_lim"])) < 10^(-10) ){
        if ( ncps_found == 0 ){
          S <- rbind(S, c(interval, rep(NA, ncol(S) - 2)))
        } else if ( ncps_found == max_cps ){
          S <- rbind(S, c(interval, b2[1:ncps_found], results_phi$results$d[1:ncps_found]))
        } else if ( ncps_found <= max_cps ){
          S <- rbind(S, c(interval, b2[1:ncps_found], rep(NA, max_cps - ncps_found), results_phi$results$d[1:ncps_found], rep(NA, max_cps - ncps_found)))
        } else {
          S <- cbind(S[, 1:(2 + max_cps), drop=FALSE], matrix(NA, nrow=nrow(S), ncol=ncps_found - max_cps), S[, -(1:(2 + max_cps)), drop=FALSE],
                 matrix(NA, nrow=nrow(S), ncol=ncps_found - max_cps))
          colnames(S) <- c("lower_lim", "upper_lim", paste0("b", 1:((ncol(S) - 2)/2)), paste0("d", 1:((ncol(S) - 2)/2)))
          S <- rbind(S, c(interval, b2[1:ncps_found], results_phi$results$d[1:ncps_found]))
          max_cps <- ncps_found
        }
        eps <- eps0
      } else {
        eps <- eps/2
      }

    }


  } else {

    # Calculate interval s.t. the given values of b & d are obtained
    interval <- calculate_interval(results, nu=nu, nu2=nu2, phi_obs=phi_obs, h1=h1, h2=h2, b0=tau)
    if ( length(b) >= 1 ){
      S <- matrix(c(interval, as.matrix(c(b, d))), nrow=1)
      colnames(S) <- c("lower_lim", "upper_lim", paste0("b", 1:length(b)), paste0("d", 1:length(d)))
    } else {
      S <- matrix(c(interval, as.matrix(c(NA, NA))), nrow=1)
      colnames(S) <- c("lower_lim", "upper_lim", "b1", "d1")
    }

    max_cps <- max(length(b), 1)

    # Find other intervals
    eps <- eps0
    upper_limit <- ifelse(model == "var", 1, Inf)
    while ( max(S[,"upper_lim"]) < upper_limit ){
      if ( model == "mean" ){
        phi <- max(S[,"upper_lim"]) + eps
        x_phi <- calculate_x_phi(x, nu, phi, nu2=nu2, phi_obs=phi_obs)
      } else if ( model == "slope" ){
        phi <- max(S[,"upper_lim"]) + eps
        x_phi <- calculate_x_phi(x, nu, phi, phi_obs=phi_obs)
      } else if ( model == "var" ){
        phi <- min(max(S[,"upper_lim"]) + eps, 1)
        x_phi <- calculate_x_phi_var(x, tau, phi, phi_obs, h1, h2)
      }
      results_phi <- find_changepoints(x_phi, method, cp_params, model)
      b2 <- results_phi$changepoints
      
      interval <- calculate_interval(results_phi, x=x, nu=nu, nu2=nu2, phi_obs=phi_obs, h1=h1, h2=h2, b0=tau)

      # Check this interval is the next one
      if ( abs(interval[1] - max(S[,"upper_lim"])) < 10^(-10) ){
        if ( is.na(b2[1]) ){
          S <- rbind(S, c(interval, rep(NA, ncol(S) - 2)))
        } else if ( length(b2) == max_cps ){
          S <- rbind(S, c(interval, b2, results_phi$results$d))
        } else if ( length(b2) <= max_cps ){
          S <- rbind(S, c(interval, b2, rep(NA, max_cps - length(b2)), results_phi$results$d, rep(NA, max_cps - length(b2))))
        } else {
          S <- cbind(S[, 1:(2 + max_cps), drop=FALSE], matrix(NA, nrow=nrow(S), ncol=length(b2) - max_cps), 
                     S[, -(1:(2 + max_cps)), drop=FALSE], matrix(NA, nrow=nrow(S), ncol=length(b2) - max_cps))
          colnames(S) <- c("lower_lim", "upper_lim", paste0("b", 1:((ncol(S) - 2)/2)), paste0("d", 1:((ncol(S) - 2)/2)))
          S <- rbind(S, c(interval, b2, results_phi$results$d))
          max_cps <- length(b2)
        }
        eps <- eps0
      } else {
        eps <- eps/2
      }
      
    }

    eps <- eps0
    lower_limit <- ifelse(model == "var", 0, -Inf)
    while ( min(S[,"lower_lim"]) > lower_limit ){
      if ( model == "mean" ){
        phi <- min(S[,"lower_lim"]) - eps
        x_phi <- calculate_x_phi(x, nu, phi, nu2=nu2, phi_obs=phi_obs)
      } else if ( model == "slope" ){
        phi <- min(S[,"lower_lim"]) - eps
        x_phi <- calculate_x_phi(x, nu, phi, phi_obs=phi_obs)
      } else if ( model == "var" ){
        phi <- max(min(S[,"lower_lim"]) - eps, 0)
        x_phi <- calculate_x_phi_var(x, tau, phi, phi_obs, h1, h2)
      }
      results_phi <- find_changepoints(x_phi, method, cp_params, model)
      b2 <- results_phi$changepoints
      
      interval <- calculate_interval(results_phi, x=x, nu=nu, nu2=nu2, phi_obs=phi_obs, h1=h1, h2=h2, b0=tau)

      # Check this interval is the next one
      if ( abs(interval[2] - min(S[,"lower_lim"])) < 10^(-10) ){
        if ( is.na(b2[1]) ){
          S <- rbind(S, c(interval, rep(NA, ncol(S) - 2)))
        } else if ( length(b2) == max_cps ){
          S <- rbind(S, c(interval, b2, results_phi$results$d))
        } else if ( length(b2) <= max_cps ){
          S <- rbind(S, c(interval, b2, rep(NA, max_cps - length(b2)), results_phi$results$d, rep(NA, max_cps - length(b2))))
        } else {
          S <- cbind(S[, 1:(2 + max_cps), drop=FALSE], matrix(NA, nrow=nrow(S), ncol=length(b2) - max_cps), S[, -(1:(2 + max_cps)), drop=FALSE],
                 matrix(NA, nrow=nrow(S), ncol=length(b2) - max_cps))
          colnames(S) <- c("lower_lim", "upper_lim", paste0("b", 1:((ncol(S) - 2)/2)), paste0("d", 1:((ncol(S) - 2)/2)))
          S <- rbind(S, c(interval, b2, results_phi$results$d))
          max_cps <- length(b2)
        }
        eps <- eps0
      } else {
        eps <- eps/2
      }

    }

  }

  S <- data.frame(S)
  colnames(S) <- c("lower_lim", "upper_lim", paste0("b", 1:max_cps), paste0("d", 1:max_cps))

  return(S)
}
