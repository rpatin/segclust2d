# neighbors
#' neighbors tests whether neighbors of point k,P can  be used to re-initialize the EM algorithm and to improve the log-likelihood.
#' @param x the initial dataset
#' @param L the likelihood
#' @param k the points of interest 
#' @param P the number of class
#' @param lmin minimal size of the segment to be implemented
#' @return smoothin likelihood
#' 
neighbors <- function (x, L,k,param,P,lmin, eps,sameSigma) {
  
  Kmax <- length(param)
  
  # first left finite neighbor
  firstLeftFiniteNeighbour <- k-1
 while(!is.finite(L[  firstLeftFiniteNeighbour]) &    firstLeftFiniteNeighbour >= P){
   firstLeftFiniteNeighbour <-   firstLeftFiniteNeighbour -1   
 }
  a  = (L[(  firstLeftFiniteNeighbour)])
  Kleft =   firstLeftFiniteNeighbour
  
  
  if ( Kleft<P ){
    phiLeft          = initialisePhi(P=P)
    out.EMLeft = list(lvinc=- Inf)
  }    else { 
    phiLeft  = param[[Kleft]]$phi     
    G        = Gmixt_simultanee(Don = x,lmin=lmin,phi=phiLeft) ## computes the cost matrix
    
    out.DP     = DynProg(G,k) ## produces the best segmentation with the given cost matrix in k segment
    t.est      = out.DP$t.est
    J.est      = out.DP$J.est
    ruptLeft      = matrix(ncol=2,c(c(1,t.est[k,1:(k-1)]+1),t.est[k,]))
    out.EMLeft    = EM.algo_simultanee(x = x,rupt = ruptLeft,P = P,phi = phiLeft, eps,sameSigma)
  } #end else
  
  
  
  # max left neighbor
  a  = max(L[1:(k-1)])
  K1 = which.max(L[1:(k-1)])
  
  
  if ( (is.null(K1)) || (a==-Inf) || (is.na(a) ) ){
    K1            =-Inf
    phi1          = initialisePhi(P=P)
    out.EM1 = list(lvinc=- Inf)
  }    else { 
    phi1                     = param[[K1]]$phi     
    G                        = Gmixt_simultanee(Don = x,lmin=lmin,phi=phi1) ## computes the cost matrix
    
    out.DP     = DynProg(G,k) ## produces the best segmentation with the given cost matrix in k segment
    t.est      = out.DP$t.est
    J.est      = out.DP$J.est
    rupt1      = matrix(ncol=2,c(c(1,t.est[k,1:(k-1)]+1),t.est[k,]))
    out.EM1    = EM.algo_simultanee(x = x,rupt = rupt1,P = P,phi = phi1, eps,sameSigma)
  } #end else
  
  
  # first right finite neighbor
  firstRightFiniteNeighbour <- k+1
  while(!is.finite(L[  firstRightFiniteNeighbour]) &    firstRightFiniteNeighbour <=Kmax){
    firstRightFiniteNeighbour <-   firstRightFiniteNeighbour +1  
  }
  a  = (L[(  firstRightFiniteNeighbour)])
  KRight =   firstRightFiniteNeighbour
  
  
  if ( KRight>Kmax ){
    phiRight          = initialisePhi(P=P)
    out.EMRight = list(lvinc=- Inf)
  }    else { 
    phiRight  = param[[KRight]]$phi     
    G        = Gmixt_simultanee(Don = x,lmin=lmin,phi=phiRight) ## computes the cost matrix
    
    out.DP     = DynProg(G,k) ## produces the best segmentation with the given cost matrix in k segment
    t.est      = out.DP$t.est
    J.est      = out.DP$J.est
    ruptRight      = matrix(ncol=2,c(c(1,t.est[k,1:(k-1)]+1),t.est[k,]))
    out.EMRight    = EM.algo_simultanee(x = x,rupt = ruptRight,P = P,phi = phiRight, eps,sameSigma)
  } #end else
  
  
  
  # right max neighbor 
  a  = max(L[(k+1):Kmax])  
  K2 = which.max(L[(k+1):Kmax])
  K2 = K2 + k
  
  if ( (K2==0) || (a==-Inf) || (is.na(a)) ){
    K2            = -Inf
    phi2          = initialisePhi(P=P)
    out.EM2       = list(lvinc=-Inf)
  }   else {
    phi2                     = param[[K2]]$phi
    G                        = Gmixt_simultanee(x,lmin=lmin,phi2)

    out.DP     = DynProg(G,k)
    t.est      = out.DP$t.est
    J.est      = out.DP$J.est
    rupt2      = matrix(ncol=2,c(c(1,t.est[k,1:(k-1)]+1),t.est[k,]))
    out.EM2    = EM.algo_simultanee(x,rupt2,P,phi2, eps,sameSigma)
    
  } #end else
  
  
  
  
  # choice of the best likelihood
  a = which.max( c(out.EM1$lvinc, out.EMLeft$lvinc, out.EMRight$lvinc, out.EM2$lvinc,  L[k]) ) 
  
  
  # parameter update
  if (a==1){
    param[[k]] = list(phi = out.EM1$phi, rupt = rupt1, tau=out.EM1$tau, cluster=apply(out.EM1$tau, 1, which.max))
    L[k] = out.EM1$lvinc}
  if (a==2){
    param[[k]] = list(phi = out.EMLeft$phi, rupt = ruptLeft, tau=out.EMLeft$tau, cluster=apply(out.EMLeft$tau, 1, which.max))
    L[k] = out.EMLeft$lvinc}
  if (a==3){
    param[[k]] = list(phi = out.EMRight$phi, rupt = ruptRight, tau=out.EMRight$tau, cluster=apply(out.EMRight$tau, 1, which.max))
    L[k] = out.EMRight$lvinc}
  if (a==4) {
    param[[k]] = list(phi = out.EM2$phi,rupt = rupt2, tau=out.EM2$tau, cluster=apply(out.EM2$tau, 1, which.max))
    L[k] = out.EM2$lvinc}
  
  invisible(list(L=L,param=param))  
    
} #end function

