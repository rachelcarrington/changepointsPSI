# Contains: get_h1_h2, get_S2

get_h1_h2 <- function(h, b, n, tau=b[1], cp_bound=TRUE, gamma=1){

  # If h = NULL, use changepoint on either side to determine h1 and h2
  if ( is.null(h) ){
    cps <- c(0, sort(b, decreasing=FALSE), n)
    j <- (1:length(cps))[cps==tau]
    h1 <- ceiling(gamma * (cps[j] - cps[j - 1]))
    h2 <- ceiling(gamma * (cps[j + 1] - cps[j]))

  # If h = (h1, h2)
  } else if ( length(h) >= 2 ) {
    if ( cp_bound ){
      cps <- c(0, sort(b[b!=tau]), n)
      if ( sum(cps %in% (tau - h[1] + 1):tau) >= 1 ){
        h1 <- tau - max(cps[cps %in% (tau - h[1] + 1):tau])
      } else {
        h1 <- h[1]
      }
      if ( sum(cps %in% (tau + 1):(tau + h[2])) >= 1 ){
        h2 <- max(cps[cps %in% (tau + 1):(tau + h[2])]) - tau
      } else {
        h2 <- h[2]
      }
    } else {
      if ( tau - h[1] < 0 ){
        h1 <- tau
        h2 <- h[2]
      } else if ( tau + h[2] > n ){
        h1 <- h[1]
        h2 <- n - tau
      } else {
        h1 <- h[1]
        h2 <- h[2]
      }
    }

  # If h is a scalar
  } else {
    stopifnot( h >= 2 )
    if ( cp_bound ){
      cps <- c(0, sort(b[b!=tau]), n)
      if ( sum(cps %in% (tau - h + 1):tau) >= 1 ){
        h1 <- tau - max(cps[cps %in% (tau - h + 1):tau])
      } else {
        h1 <- h
      }
      if ( sum(cps %in% (tau + 1):(tau + h)) >= 1 ){
        h2 <- max(cps[cps %in% (tau + 1):(tau + h)]) - tau
      } else {
        h2 <- h
      }
    } else {
      if ( tau - h < 0 ){
        h1 <- tau
      } else {
        h1 <- h
      }
      if ( tau + h > n ){
        h2 <- n - tau
      } else {
        h2 <- h
      }
    }
  }

  return(c(h1, h2))
}



get_S <- function(S_all, h, b, tau=b[1]){
  # Calculate S from S_all
  
  # h is the window size; b is the vector of detected changepoints; tau is the changepoint of interest

  max_cps_found <- (ncol(S_all) - 2)/2
  inS <- rep(1, nrow(S_all))

  # Option 1. If h = NULL, S should contain intervals which contain the correct combination of changepoints
  if ( is.null(h) ){
    k <- length(b)
    # First, remove intervals which contain too many or too few changepoints
    if ( max_cps_found > k ){
      inS[!is.na(S_all[,(k + 3)])] <- 0 # remove rows where the (k+1)th changepoint is not NA (too many)
    }
    inS[is.na(S_all[,(k + 2)])] <- 0 # remove rows where the (k)th changepoint is NA (too few)
    
    # Now, for remaining intervals (which have the correct # of CPs), check they are the correct ones
    cps <- sort(b)
    if ( k == 1 ){
      inS[S_all[,3] != b[1]] <- 0
    } else {
      inS[inS == 1][colSums((apply(S_all[inS == 1, 3:(k + 2), drop=FALSE], 1, sort)) - cps) != 0] <- 0
    }

    S <- S_all[inS == 1,,drop=FALSE]
    
  # Option 2. If h is not NULL, S2 contains intervals which contain the changepoint of interest
  } else {
    S <- S_all[S_all$b1 == tau, ]
    if ( max_cps_found > 1 ){
      for ( j in 2:max_cps_found ){
        S <- rbind(S, S_all[S_all[,j + 2] == tau, ])
      }
    }
    S <- S[!is.na(S$b1),,drop=FALSE]
  }
  
  return(S)
}