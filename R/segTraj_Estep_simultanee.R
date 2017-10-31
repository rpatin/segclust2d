# Estep_simultanee
#' Estep_simultanee computes posterior probabilities and incomplete-data  log-likelihood for mixture models
#' @param logdensity is a  K*P matrix containing the conditinal log-densities for each segment
#' @param phi   a list containing the parameters of the mixture
#' @param eps eps
#' @return a list with tau a K*P matrix, tau kj is the posterior probability for segment k to belong to classe j and lvinc, the incomplete log vrais P(X=x)

Estep_simultanee <- function (logdensity,phi, eps=1e-9){
  K = nrow(logdensity)
  P = ncol(logdensity)
  tau     = matrix((log( phi$prop)),K,P,byrow=TRUE)+logdensity
  ## avoiding empty clusters, adding pseudo counts
  # listeEmpty = which(colSums(tau) <eps)
  # if(length(listeEmpty)>1){
  #   tau[,listeEmpty] <- tau[,listeEmpty]+eps
  #   tau <- sweep(tau,MARGIN = 1, STATS = rowSums(tau), FUN = '/')
  # }
  tau_max = apply(tau,1,max)
  tau     = exp(tau-matrix(tau_max,K,P))
  tau     = sweep(tau,1,STATS = apply(tau,1,sum), FUN="/")

  lvinc   = sum(log( apply(tau,1,sum)) + tau_max) #sum(log P(Xk=xk))
  invisible(list(tau,lvinc))
}

