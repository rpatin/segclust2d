# EM.algo_simultanee
#' EM.algo_simultanee caculates the MLE of phi for given change-point instants
# and for a fixed number of clusters
#' @param rupt ahe sequence of change points
#' @param P   number of clusters
#' @param  phi starting value for the  parameter
#' @return a list with  phi, the MLE, tau =(taukj) the probability for segment k to belong to classe,lvinc = lvinc,empty = empty,dv = dv
#' @export
 
EM.algo_simultanee <- function(x,rupt,P,phi, eps=1e-6,sameSigma=FALSE){

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
    logdensity = t( apply(rupt,1,FUN=function(y) logdens_simultanee(   x[, y[1]:y[2] ],phi)))

    Estepout   = Estep_simultanee(logdensity,phi)
    tau        = Estepout[[1]]
    
    
    lvinc      = Estepout[[2]]

    phi        = Mstep_simultanee(x,rupt,tau,phi,sameSigma)
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
