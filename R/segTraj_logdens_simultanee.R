# logdens_simultanee
#' logdens_simultanee calculates Gaussian log-densities for a bivariate signal under a mixture model with P components
#' @param xk  the bivaraiate signal
#' @param phi   parameters of the mixture, P components
#' @return the value of the log density
logdens_simultanee <- function(xk,phi){
  
  P <- length(phi$prop)
  m  = phi$mu
  s  = phi$sigma
  
  nk = dim(xk)[2]
  tmp = matrix(ncol=P,nrow=1)
  
  tmp=sapply(1:P, function(p){
    -nk*log(sqrt(2*pi)*s[1,p])- 0.5*sum (    (xk[1,]  - m[1,p])^2  )/s[1,p]^2-nk*log(sqrt(2*pi)*s[2,p])- 0.5*sum (    (xk[2,]  - m[2,p])^2  )/s[2,p]^2
  } )
  
  invisible(tmp)
  
}
