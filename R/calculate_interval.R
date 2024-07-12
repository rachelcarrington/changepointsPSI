#' Calculate interval of \eqn{\phi} values which return a given set of changepoints.
#'
#' @description For a given changepoint algorithm, model, and dataset, calculate the interval of \eqn{\phi} values such that applying the algorithm to \eqn{X'(\phi)} will give 
#' a particular set of changepoint estimates and directions of change. (Wrapper function for \code{calculate_interval_bs}, \code{calculate_interval_wbs}, and \code{calculate_interval_not}.)
#'
#' @param results Output of \code{find_changepoints}.
#' @param x Vector of data (if different from that supplied to \code{find_changepoints}).
#' @param b0 Changepoint of interest; only needed if \code{model = "var"}, as otherwise this information is captured by the vector \code{nu}.
#' If this is not supplied, the first detected changepoint will be used.
#' @param nu Vector which defines the test statistic; not required if \code{model == "var"}. See below.
#' @param nu2 Value of \eqn{||\nu||_2^2} (optional).
#' @param phi_obs Observed value of the test statistic (optional).
#' @param h1 Left window size; only needed if \code{model = "var"}.
#' @param h2 Right window size; only needed if \code{model = "var"}.
#' @param maxiter Maximum number of changepoints to detect.
#' @param C0 Vector containing cumulative sums of \code{x} (if \code{model = "mean"} or code{model = "slope"}) or \code{x^2} (if \code{model = "var"}).
#' @param C0_nu Vector containing cumulative sums of \code{nu}.
#' @param autocor Logical: \code{TRUE} if data is autocorrelated, otherwise \code{FALSE} (default).
#' @param Sigma Covariance matrix of the data (only needed if \code{autocor = TRUE}).
#'
#' @return A 2-dimensional vector.
#' @export
#'
#' @details
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
#' The window sizes \code{h1} and \code{h2} must be supplied to the function.
#'
#' The parameters \code{C0}, \code{nu2}, and \code{C0_nu} are optional, but will speed up computation time if the function is used multiple times (e.g. within \code{calculate_pvals_all}).
#' The latter two are only relevant if \code{model = "mean"} or \code{model = "slope"}.
#' 
#' @examples
#' set.seed(100)
#' x <- rnorm(100) + c(rep(1,50), rep(-1,50))
#' results <- find_changepoints(x, "bs", params=list(threshold=4))
#' b <- results$changepoints
#' h <- 10
#' nu <- c(rep(0, b[1] - h), rep(1/h, h), rep(-1/h, h), rep(0, length(x) - b[1] - h))
#' calculate_interval(results, nu=nu)
#'
calculate_interval <- function(results, x=NULL, b0=NULL, nu=NULL, nu2=NULL, phi_obs=NULL, h1=NULL, h2=NULL, maxiter=NULL, C0=NULL, C0_nu=NULL, autocor=FALSE, Sigma=NULL ){
  
  if ( is.null(x) ){
    x <- results$x
  }
  method <- results$params$method
  model <- results$params$model
  b <- results$changepoints
  d <- results$results$d
  threshold <- results$params$threshold
  if ( is.null(maxiter) ){
    maxiter <- results$params$maxiter
  }
  
  # For change in variance model, if loss = "icss" then use ICSS loss function
  icss <- FALSE
  if ( length(results$params$loss) == 1 ){
    if ( results$params$loss == "icss" ){
      icss <- TRUE
    }
  }

  if ( method == "bs" ){

    interval <- calculate_interval_bs(x, b, d, b0=b0, nu=nu, threshold=threshold, phi_var=nu2, phi_obs=phi_obs, maxiter=maxiter, model=model, h1=h1, h2=h2, C0=C0, 
                                      C0_nu=C0_nu, icss=icss, autocor=autocor, Sigma=Sigma)

  } else if ( method == "wbs" ){
    
    s <- results$results$s
    e <- results$results$e
    random_samples <- results$params$random_samples
    interval <- calculate_interval_wbs(x, b=b, d=d, s=s, e=e, random_samples=random_samples, nu=nu, phi_var=nu2, phi_obs=phi_obs, model=model, maxiter=maxiter, 
                                       threshold=threshold, tau0=b0, h1=h1, h2=h2, icss=icss, autocor=autocor, Sigma=Sigma)
    
  } else if ( method=="not" ){

    s <- results$results$s
    e <- results$results$e
    random_samples <- results$params$random_samples
    interval <- calculate_interval_not(x, b=b, d=d, s=s, e=e, random_samples=random_samples, nu=nu, phi_var=nu2, phi_obs=phi_obs, model=model, maxiter=maxiter, 
                                       threshold=threshold, tau0=b0, h1=h1, h2=h2, icss=icss, autocor=autocor, Sigma=Sigma)
  
  }

  return(interval)

}