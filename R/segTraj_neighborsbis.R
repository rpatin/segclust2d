# neighbors
#' neighbors tests whether neighbors of point k,P can  be used to re-initialize the EM algorithm and to improve the log-likelihood.
#' @param x the initial dataset
#' @param L the likelihood
#' @param k the points of interest
#' @param P the number of class
#' @param lmin minimal size of the segment to be implemented
#' @param kv.hull convex hull of likelihood
#' @param param param outputs of segmentation
#' @param eps eps
#' @param sameSigma should segments have same variance ?
#' @param pureR should algorithm use only R functions or benefit from Rcpp faster algorithm
#' @return smoothin likelihood
#'
neighborsbis <- function (kv.hull,x, L,k,param,P,lmin, eps,sameSigma=TRUE, pureR = F) {

  for (j in 1:length(kv.hull)){
    K1=kv.hull[j]
    a  = L[K1]
    if (a==-Inf){
      K1            =-Inf
      phi1          = initialisePhi(P=P)
      out.EM1 = list(lvinc=- Inf)
    }    else {
      phi1                     = param[[K1]]$phi
      if(pureR){
        G = Gmixt_simultanee(x,lmin,phi1) ## computes the cost matrix
        out.DP     = DynProg(G,k) ## produces the best segmentation with the given cost matrix in k segment
      } else {
        G = Gmixt_simultanee_fullcpp(x, lmin=lmin, phi1$prop, phi1$mu, phi1$sigma)
        out.DP     = wrap_dynprog_cpp(G,k)
      }
      t.est      = out.DP$t.est
      J.est      = out.DP$J.est
      rupt1      = matrix(ncol=2,c(c(1,t.est[k,1:(k-1)]+1),t.est[k,]))
      if(pureR){
        out.EM1 = EM.algo_simultanee(x = x,rupt = rupt1,P = P,phi = phi1, eps,sameSigma)
      } else {
        out.EM1 = EM.algo_simultanee_Cpp(x = x,rupt = rupt1,P = P,phi = phi1, eps,sameSigma)
      }
    } #end else
   if (out.EM1$lvinc>L[k])
     {
     param[[k]] = list(phi = out.EM1$phi, rupt = rupt1, tau=out.EM1$tau, cluster=apply(out.EM1$tau, 1, which.max))
     L[k] = out.EM1$lvinc
   }
   }

  invisible(list(L=L,param=param))

} #end function

