# initialisePhi
#' initialisePhi is the constructor for a set of parameters for a segclust model
#' P number of classes
#' val the value ti be used for tinitiliastation default is -Inf
#' @return a set of paraemter phi
initialisePhi <- function(P, val=-Inf)
{
  return(list(mu=matrix(val, nrow = 2, ncol=P), sigma=matrix(val, nrow = 2, ncol=P ), prop=rep(1/P, P) ))
}