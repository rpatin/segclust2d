# EM.algo_simultanee C++
#' EM.algo_simultanee caculates the MLE of phi for given change-point instants
#'  and for a fixed number of clusters
#' @param rupt ahe sequence of change points
#' @param P   number of clusters
#' @param phi starting value for the  parameter
#' @param x bivariate signal
#' @param eps eps
#' @param sameSigma TRUE if segments have the same variance
#' @return a list with  phi, the MLE, tau =(taukj) the probability for segment k to belong to classe,lvinc = lvinc,empty = empty,dv = dv
#' @export

EM.algo_simultanee_Cpp <- function(x,rupt,P,phi, eps=1e-6,sameSigma=FALSE){

  K     = nrow(rupt)
  delta = 1
  empty = 0
  dv    = 0
  tau   = matrix(1,nrow = K,ncol = P)
  iter  = 0
  np    = apply(tau,2,sum)

  while ( (delta>=1e-4) & (min(np)>eps) & (iter<=500 )){

    iter       = iter+1
    phi_temp   = phi
    # logdensity = t( apply(rupt,1,FUN=function(y) logdens_simultanee(   x[, y[1]:y[2] ],phi)))
    logdensity = t( apply(rupt,1,FUN=function(y) logdens_simultanee_cpp(   x[, y[1]:y[2] ], phi$mu, phi$sigma, phi$prop)))

    Estepout   = Estep_simultanee(logdensity,phi)
    tau        = Estepout[[1]]


    lvinc      = Estepout[[2]]

    phi        = Mstep_simultanee_cpp(x,rupt,tau,phi,sameSigma)
    np         = apply(tau,2,sum)

    delta      =max(unlist(lapply(names(phi),function(d) {max(abs(phi_temp[[d]]-phi[[d]])/phi[[d]])})))
  }

  if (min(np)<eps){
    empty = 1
    lvinc = -Inf
  }

  if (iter>5000){
    dv    = 2
    lvinc = -Inf
  }

  rm(delta,logdensity)


  invisible(list(phi = phi,tau = tau,lvinc = lvinc,empty = empty,dv = dv))

}


# Mstep_simultanee C++
#' Mstep_simultanee computes the MLE within the EM framework

#' @param x the bivariate signal
#' @param rupt the rupture dataframe
#' @param phi the parameters of the mixture
#' @param tau the K*P matrix containing posterior probabilities of membership to clusters
#' @param sameSigma whether segments have the same variance
#' @return phi the updated value of the pramaeters

Mstep_simultanee_cpp <- function (x,rupt,tau,phi,sameSigma=TRUE) {

  K = nrow(tau)
  P = ncol(tau)
  m = matrix(nrow=2,ncol=P)
  s = matrix(nrow=2,ncol=P)
  prop = matrix(nrow=1,ncol=P)
  # Yk = apply(rupt,1,FUN=function(y) rowSums(x[,y[1]:y[2]]))
  Yk = apply_rowSums(rupt,x)
  rownames(Yk) <- rownames(x)

  nk = rupt[,2]-rupt[,1]+1
  n  = sum(nk)

  #
  np    = nk %*% tau
  m=Yk %*% tau/rep(np,each=2)
  if(!sameSigma){
    for (i in 1:2){
      # s[i,]=  colSums( tau*(sapply(1:P, function(p) {apply(rupt,1,FUN=function(y) sum((x[i,y[1]:y[2]]-m[i,p])^2   ))})))
      s[i,]  <-  colsums_sapply(i, rupt, x, m, tau)
    }
    s=sqrt(s/rep(np,each = 2))
  } else
  {
    for (i in 1:2){
      s[i,]=  rep(sum( tau*(sapply(1:P, function(p) {apply(rupt,1,FUN=function(y) sum((x[i,y[1]:y[2]]-m[i,p])^2   ))}))), P)
    }
    s=sqrt(s/n)

  }

  # prop = apply(tau,2,sum)/K
  # emptyCluster = which(prop==0)
  # if(length(emptyCluster)>0){
  #   prop = pmax(prop, eps)
  #   prop = prop /sum(prop)
  #   for (d in emptyCluster){
  #     m[,d]=rep(0,2)
  #     s[,d]=rep(1e9,2)
  #   }
  # }

  prop = apply(tau,2,sum)/K
  b    = order(m[1,])
  m    = m[,b]
  s    = s[,b]
  prop = prop[b]
  phi  = list(mu=m,sigma=s,prop=prop)

  invisible(phi)
}
