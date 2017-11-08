# initialisePhi
#' initialisePhi is the constructor for a set of parameters for a segclust model
#' @param P number of classes
#' @param val the value used for initilisation default is -Inf
#' @return a set of parameter phi
initialisePhi <- function(P, val=-Inf)
{
  return(list(mu=matrix(val, nrow = 2, ncol=P), sigma=matrix(val, nrow = 2, ncol=P ), prop=rep(1/P, P) ))
}
