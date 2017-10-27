# DynProg RCpp
#' DynProg computes the change points given a cost matrix matD
#' and a maximum number of segments Kmax
#' @param matD the cost Matrix os size n x n
#' @param Kmax the maximal number of segments
#' @return a list with J.est a vector with Kmax value, the Kth is the minimum contrast for a model with K segments (-J.est is the log-likelihood) and with  t.test a matrix, line K are the coordinates of the change points for a model with K segments
#' @export
#'
wrap_dynprog_cpp <- function(G,K){
  out <- DynProg_algo_cpp(G,K)
  out$J.est = out$J.est[,1]
  out$t.est[which(is.na(out$t.est))]=Inf

  return(out)
}


# DynProg
#' DynProg computes the change points given a cost matrix matD
#' and a maximum number of segments Kmax
#' @param matD the cost Matrix os size n x n
#' @param Kmax the maximal number of segments
#' @return a list with J.est a vector with Kmax value, the Kth is the minimum contrast for a model with K segments (-J.est is the log-likelihood) and with  t.test a matrix, line K are the coordinates of the change points for a model with K segments
#' @export
DynProg<-function(matD,Kmax){

  N<-dim(matD)[1]
  if (Kmax>N){
    cat("Kmax ", Kmax, "is greater than N ",N,"\n")
    cat("Kmax is supposed to be equal to N :", N,"\n")
  }

  I<-matrix(Inf,Kmax,N)
  t<-matrix(0,Kmax,N)
  I[1,]=matD[1,]
  matD=t(matD) ## a quoi ca sert ?

  if (Kmax>2) {

    for (k in 2:(Kmax-1)){
      for (L in k:N){
        I[k,L]<-min(I[(k-1),1:(L-1)]+matD[L,2:L])

        if(I[k,L]!=Inf){
          t[(k-1),L]<-which.min(I[(k-1),1:(L-1)]+matD[L,2:L])
        } else {
          t[(k-1),L] = Inf
        } # end else
      } #end L
    }#end k
  } #end K


  I[Kmax,N]<-min(I[(Kmax-1),1:(N-1)]+matD[N,2:N])

  if(I[Kmax,N]!=Inf){
    t[(Kmax-1),N]<-which.min(I[(Kmax-1),1:(N-1)]+matD[N,2:N])
  } #end if


# *** computation of breakpoint instants ***

  t.est<-matrix(0,Kmax,Kmax)
  diag(t.est)<-N
  for (K in 2:Kmax){
    for (k in seq(K-1,1,by=-1)){
      if(t.est[K,k+1]!=0){
        t.est[K,k]<-t[k,t.est[K,k+1]]
      } #end if
    } #end k
  } #end K
  t.est[which(is.na(t.est))]=Inf
  list(J.est = I[,N],t.est = t.est)
  #cat("Kmax", Kmax,"\n")
}

