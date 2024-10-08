% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/l0_segmentation_psi.R
\name{l0_segmentation_psi}
\alias{l0_segmentation_psi}
\title{Post-selection inference for the change in mean model, using L0 segmentation to estimate changepoints.}
\usage{
l0_segmentation_psi(
  x,
  lambda,
  h,
  N = 1,
  sigma2 = NULL,
  sig = 4,
  include_original = TRUE,
  num_pvals = NULL,
  cp_bound = TRUE,
  gamma = 1
)
}
\arguments{
\item{x}{Numeric vector of data.}

\item{lambda}{Changepoint detection threshold.}

\item{h}{Window size.}

\item{N}{Number of \eqn{\psi} samples to take (defaults to 1).}

\item{sigma2}{Variance of \code{x}.}

\item{sig}{Tuning parameter.}

\item{include_original}{Logical; whether to include observed value as \eqn{\psi} in place as one of the random samples;
defaults to \code{TRUE}.}

\item{num_pvals}{Maximum number of p-values to calculate.}
}
\value{
A list:
\itemize{
\item \code{b} Vector of changepoints
\item \code{p_value} Vector of p-values
\item \code{p_value_orig} Vector of p-values obtained using fixed \eqn{\psi = \psi_{obs}}
\item \code{phi_obs} The observed value of the test statistic.
\item \code{P_both} Matrix containing values of \eqn{Pr(|\phi| > |\phi_{obs}| \& \phi \in S)}
\item \code{P_phi_in_S} Matrix containing values of \eqn{Pr(\phi \in S)}
\item \code{P_both_orig} Vector containing values of \eqn{Pr(|\phi| > |\phi_{obs}| \& \phi \in S)} for fixed \eqn{\psi = \psi_{obs}}
\item \code{P_phi_in_S_orig} Vector containing values of \eqn{Pr(\phi \in S)} for fixed \eqn{\psi = \psi_{obs}}
}
}
\description{
Post-selection inference for L0 segmentation. This uses the functions \code{changepoint_estimates} and
\code{changepoint_inference} from the package \code{ChangepointInference}.
}
\examples{
set.seed(100)
x <- rnorm(100) + c(rep(1,40), rep(-1,20), rep(1,40))
l0_segmentation_psi(x, 4, 10, 10)

}
