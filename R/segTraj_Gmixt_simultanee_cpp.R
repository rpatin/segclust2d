# Gmixt_simultanee
#'
#' Gmixt_simultanee calculates the cost matrix for a segmentation/clustering model
#' @param Don the bivariate  signal
#' @param lmin the minimum size for a segment
#' @param phi  the  parameters of the mixture
#' @return a matrix  G(i,j), the mixture density for segment between points (i+1) to j
#'         = deqn{\sum_{p=1}^P \log (\pi_p f(y^{ij};\theta_p))}
#'          Rq: this density if factorized in order to avoid numerical zeros in the log

Gmixt_simultanee_cpp <- function(Don,lmin,phi){

  P = length(phi$prop)
  m    = phi$mu
  s    = phi$sigma
  prop = phi$prop

  n = dim(Don)[2]
  G = list()

  for (signal in 1:2){

    z = Don[signal,]
    lg  = lmin:n  # possible position for the end of  first segment
    zi  = cumsum(z) # z cumul
    zi=zi[lg]
    z2  = z^2
    z2i = cumsum(z2)
    z2i=z2i[lg]

    #rappel: on fait du plus court chemin
    #        donc on prend -LV
    # G[[signal]][1, l] contains the likelihood of the segment starting in 1 and finishing

    G[[signal]] <- matrix(Inf,ncol=n,nrow=n)

    ## The following code makes use of vectoriel facilities
    wk  <- (z2i/lg-(zi/lg)^2)
    #  old # dkp <- (sweep(repmat(t(zi/lg),P,1), MARGIN = 1, STATS = m[signal,]))^2
    ##  very old # Aold    = (wkold+dkpold)/repmat(s[signal,]^2,1,n-lmin+1)+ log(2*pi*repmat(s[signal,]^2,1,n-lmin+1))
    tmp_repmat <- arma_repmat_transpose_divide(zi,lg,P,1)
    dkp <- (sweep_col_plus(tmp_repmat,-m[signal,]))^2

    # A <- sweep(dkp, MARGIN = 2, STATS = wk, FUN = '+')
    A <- sweep_row_plus(dkp, wk)

    # A <- sweep(A, MARGIN = 1, STATS = s[signal,]^2, FUN = '/')
    A <- sweep_col_divide(A, s[signal,]^2)

    # A <- sweep(A, MARGIN = 1, STATS = log(2*pi*s[signal,]^2), FUN = '+')
    A <- sweep_col_plus(A, log(2*pi*s[signal,]^2))

    # A <- -0.5*sweep(A, MARGIN = 2, STATS = lg, FUN = '*')
    A <- -0.5*sweep_row_times(A, lg)

    # A <- sweep(A, MARGIN = 1, STATS = log(prop), FUN='+')
    A <- sweep_col_plus(A,log(prop))

    ## finally A[k,p] contains -0.5 *\sum_{l=1^(k)} (Z_l -\mu^p)^2 /(sigma_p^2) +log(pi_p)
    ## that is the density of observing segment 1:k and beeing in cluster p

    ## the normalization using A_max is used to avoid numerical issues with the exponential
    # A_max = apply(A,2,max)
    A_max = apply_col_max(A)
    # Aprov = exp(sweep(A, MARGIN = 2, STATS = A_max, FUN = '-'))
    Aprov = exp(sweep_row_plus(A,-A_max))
    G[[signal]][1,lmin:n] = -log(apply_col_sum(Aprov)) - A_max
    # return(G)
    for (i in (2:(n-lmin+1))) {
      ni  = n-i-lmin+3
      z2i = z2i[2:ni]-z2[i-1]
      zi  = zi[2:ni]-z[i-1]
      lgi = lmin:(n-i+1)
      wk = (z2i)/(lgi)-(zi/(lgi))^2
      # dkp <- (sweep(repmat(t(zi/lgi),P,1), MARGIN = 1, STATS = m[signal,]))^2

      G[[signal]][i,(i+lmin-1):n] =  Gmixt_algo_cpp(zi, lgi, P, m[signal,], wk, s[signal,], prop)
      ## problem if i=n-lmin+1
      ## Aprov[,(i+lmin-1):n] is a vector, and apply can't be used
      ## this case is postponed at the end of the loop
    }
    # i <- (n-lmin+1)
    #  Aprov <- sweep(A, MARGIN = 1, STATS = A[,i], FUN = '-')
    #  Aprov_max = apply(Aprov,2,max)
    #  Aprov     = exp(sweep(Aprov, MARGIN = 2, STATS = Aprov_max, FUN = '-'))
    #  ## problem if i=n-lmin+1
    #  ## Aprov[,(i+lmin-1):n] is a vector, and apply can't be used
    #  ## this case is postponed at the end of the loop
    #  G[[signal]][i,(i+lmin-1):n] = -log(sum(Aprov[,(i+lmin-1):n])) - Aprov_max[(i+lmin-1):n]

  }
  # return(G)
  res <- Reduce('+', G )
  invisible(res)
}
