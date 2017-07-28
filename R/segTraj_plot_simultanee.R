# bisig_plot
#' bisig_plot draws the plots of the bivaraite signal on the same plot (scale free)
#' @param x the signal to be plotted
#' @param rupt optionnal, if given add vertical lines at change points (rupt should a vector)
#' @param mu optionnal the mean of each class of segment, 
#' @param pop optionnal the cluster to whom each segment belongs to,
#' 
#' @export
#' @return no value

bisig_plot<- function(x, rupt=NULL, mu=NULL, pop=NULL, merge.seg=FALSE){
  n <- dim(x)[2]
  if(!is.null(rupt)){
    if(!is.matrix(rupt))
      rupt       = ruptAsMat(rupt)
    K <- dim(rupt)[1]
  }
  m <- rowMeans(x)
  s <- apply(x,1,sd)
  
  x.norm <- sweep(x = x, MARGIN = 1, STATS = m)
  x.norm <- sweep(x = x.norm, MARGIN = 1, STATS = s, FUN = '/')
  plot(1:n,x.norm[1,], pch=19, col=2, ylab='', xlab = '', xaxt='n', yaxt='n')
  points(1:n, x.norm[2,], pch=19, col=3)
  if(!is.null(rupt)&!(merge.seg)){
    invisible( lapply(1:K, function(d){ abline(v=rupt[d,1], lwd=1.5)}))
  }
  if(!is.null(rupt)&(merge.seg)){
    change <- which(diff(pop)!=0)
    p <- length(change)
    ruptProv <- matrix(NA, nrow=p+1, ncol=2)
    ruptProv[1, 1] <- rupt[1,1]
    ruptProv[1, 2] <- rupt[change[1],2]
    ruptProv[p+1, 1] <- rupt[change[p]+1,1]
    ruptProv[p+1, 2] <- rupt[K,2]
    
    if(p>1){
      ruptProv[2:p,] <- t(sapply(1:(p-1), function(p_){c(rupt[change[p_]+1,1], 
                                                       rupt[change[p_+1],2])}))
    }
    invisible( lapply(1:(p+1), function(d){ abline(v=ruptProv[d,1], lwd=1.5)}))
    
  }
  if(!is.null(mu) & !is.null(pop)){
    invisible( lapply(1:K, function(d){ 
      lines(c(rupt[d,1],rupt[d,2]+1), (rep(mu[1,pop[d]],2)-m[1])/s[1], col=2)
      lines(c(rupt[d,1],rupt[d,2]+1), (rep(mu[2,pop[d]],2)-m[2])/s[2], col=3)
    }))
  }
}



# matrixRupt
#' matrixRupt transforms a vector of change point into a data.frame with start and end of every segment
#' @param x the 
#' @param vectorRupt the vector containing the change point
#' 
#' @export
#' @return the matrix of change point

matrixRupt <- function(x,  vectorRupt){
 n <- ncol(x)
  posZeros <-  which(vectorRupt==0)
 if(length(posZeros)>0){
  K <- posZeros[1]-1
  vectorRupt <- vectorRupt[1:K]
 } else { }
  return(matrix(c(1, vectorRupt[1:(K-1)]+1, vectorRupt[1:(K)]), ncol=2))
}