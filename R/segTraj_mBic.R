# mBic
#' mBic computes the modified Bic criterion

#' @param Don the bivariate signal
#' @param resHSimultanee the results of hybrid_simultanee
#' @return the mBic criterion for each number of change points K

mBic <- function(Don,  resHSimultanee){
  n <- dim(Don)[2]
  I <- dim(Don)[1]
  logLength <- unlist(lapply(resHSimultanee$param, function(d){
    if(!is.null(d))
    return(sum(log(d$rupt[,2]-d$rupt[,1]+1)))
    else
      return(-Inf)
  }))
  N= I*n
 mBic <- sapply(1:length(logLength), function(k){
    if( logLength[k] == -Inf)
      return(Inf)
    else
      ((N-k*I)/2+1)* (2/N*resHSimultanee$Linc[k]+1+log(2*pi)-log(N)) 
    - 0.5*logLength[k]-(k-1)*log(n)+lgamma((N-k*I)/2+1)
})
  return (mBic)
}