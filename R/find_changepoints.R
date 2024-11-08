#' Estimate changepoints in mean, variance, or slope.
#'
#' @description Function to estimate changes in mean, variance or slope for univariate data. Wrapper function for \code{binary_segmentation}, \code{wild_binary_segmentation}, 
#' and \code{narrowest_over_threshold}. Also includes PELT using the \code{changepoint} package.
#'
#' @param x Vector of data.
#' @param method Changepoint method: one of \code{"bs"} (binary segmentation), \code{"wbs"} (wild binary segmentation), \code{"not"} (narrowest over threshold),
#' or \code{"pelt"} (PELT). For PELT the \code{changepoint} package must be installed.
#' @param params A list containing parameters for the changepoint algorithm (other than \code{x} and \code{model}). (See below.)
#' To use default values for all parameters, enter an empty list: \code{list()}.
#' @param model Changepoint model: one of \code{"mean"} (piecewise constant mean, default), \code{"slope"} (piecewise linear mean), or \code{"var"} (piecewise constant variance).
#'
#' @return A list:
#' \itemize{
#' \item \code{x} The data.
#' \item \code{results} Dataframe containing results.
#' \item \code{results_full} For wild binary segmentation or narrowest over threshold, a dataframe containing the changepoint candidates for every random interval.
#' \item \code{changepoints} Vector of changepoints detected.
#' \item \code{params} The parameters supplied to the function.
#' }
#'
#' @details 
#' The parameter \code{params} should be a list containing relevant parameters for the changepoint algorithm. This can be left empty, in which case parameters will take their
#' default values. There is no need to specify \code{x} or \code{model} as these are supplied directly to the \code{find_changepoints} function. See documentation for individual
#' functions for a list of possible parameters.
#'
#'
#' @examples
#' set.seed(100)
#' x <- rnorm(100) + c(rep(1,45), rep(-1,10), rep(1,45))
#' results_bs <- find_changepoints(x, "bs", list(threshold=3))
#' print(results_bs)
#'
#' results_wbs <- find_changepoints(x, "wbs", list(threshold=4, num_rand_samples=100))
#' print(results_wbs)
#'
#' results_not <- find_changepoints(x, "not", list(threshold=4, num_rand_samples=100))
#' print(results_not)
#'
#' results_pelt <- find_changepoints(x, "pelt", list())
#' print(results_pelt)
#'
#' x <- c(rnorm(100), rnorm(100, sd=1.5))
#' find_changepoints(x, "bs", list(threshold=4), "var")
#'
#' @export
#'
find_changepoints <- function(x, method, params=list(), model="mean"){

  stopifnot( method == "bs" || method == "wbs" || method == "not" || method == "pelt" )
  stopifnot( model == "mean" || model == "var" || model == "slope" )
  
  results_full <- NA

  if ( method == "bs" ){

    results <- do.call(binary_segmentation, c(params, model=model, x=list(x)))
    params$threshold <- results$threshold
    params$maxiter <- results$maxiter
    changepoints <- results$changepoints
    results <- results$results[results$results$cp == 1, ]

  } else if ( method == "wbs" ){

    results <- do.call(wild_binary_segmentation, c(params, model=model, x=list(x)))
    params$threshold <- results$threshold
    params$maxiter <- results$maxiter
    params$random_samples <- results$random_samples
    results_full <- results$results_full
    changepoints <- results$changepoints
    results <- results$results

  } else if ( method == "not" ){

    results <- do.call(narrowest_over_threshold, c(params, model=model, x=list(x)))
    params$threshold <- results$threshold
    params$maxiter <- results$maxiter
    params$random_samples <- results$random_samples
    results_full <- results$results_full
    changepoints <- results$changepoints
    results <- results$results

  } else if ( method == "pelt" ){ # could add in also other methods using cpt.mean or cpt.var

    if ( model == "var" ){
      results <- do.call(changepoint::cpt.var, c(params, data=list(x), method="PELT"))
      changepoints <- results@cpts[-length(results@cpts)]
      var_ests <- results@param.est$variance
      d <- ifelse(var_ests[-1] - var_ests[-length(var_ests)] > 0, 1, -1)
    } else if ( model == "mean" ){
      results <- do.call(changepoint::cpt.mean, c(params, data=list(x), method="PELT"))
      changepoints <- results@cpts[-length(results@cpts)]
      mean_ests <- results@param.est$mean
      d <- ifelse(mean_ests[-1] - mean_ests[-length(mean_ests)] > 0, 1, -1)
    }
    lrs <- s <- e <- NA
    results <- data.frame(b=changepoints, d=d)

  }
  
  return(list(x=x, results=results, results_full=results_full, changepoints=changepoints, params=c(list(model=model, method=method), params)))
}
