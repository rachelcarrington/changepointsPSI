# Estimate variance of a dataset using median absolute deviation, where data is autocorrelated.

mad_ar <- function(x, ord=1){
  
  difs <- x[-(1:ord)] - x[1:(length(x) - ord)]
  mad_est <- median(abs(difs - median(difs))) / qnorm(0.75)
  mad_est <- mad_est^2

  return(mad_est)
}