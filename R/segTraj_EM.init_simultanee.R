# EM.init_simultanee 
#' EM.init_simultanee proposes an inial value for the EM algorithm  based on a hierarchical clustering  algorithm (ascending)
#' @param x the bivariate signal
#' @param rupt the  change point instants, data.frame
# K number of segments
# P numer of clusters
# phi0 : candidate for the EM algorithm 
EM.init_simultanee <- function (x,rupt,K,P) {

m = apply(rupt,1,FUN=function(y) rowMeans(x[,y[1]:y[2]]))
v =  apply(rupt,1,FUN=function(y) apply(x[,y[1]:y[2]], 1, var))
n = apply(rupt,1,diff)+1

Dist = matrix(Inf,K,K)

for (k in (1: (K-1))){

   for (r in ( (k+1) :K)) {
   
      ybar      = (n[k]*m[,k]+n[r]*m[,r])/(n[k]+n[r])
      varpool   = sum(  n[r]*v[,r] + n[k]*v[,k] + n[r]*(m[,r]-ybar)^2 + n[k]*(m[,k]-ybar)^2 ) / (n[r]+n[k])
      Dist[k,r]    = sum(-n[k]*log(v[,k])-n[r]*log(v[,r]))+ (n[k]+n[r])*log(varpool)
      }
}
#Dist=matrice K*K = distance (associee a la log-vraisemblance) entre les segments

LR = NULL

if (nrow(Dist)>P) {
for (indice in (1:K)){
#while (length(Dist) > P){
   
   if (nrow(as.matrix(Dist)) <= P) {
   break 
   }
   
   out = which(Dist==min(Dist),arr.ind=TRUE)
   imin = out[1]
   jmin = out[2]
   Dmin = min(Dist)
   
   LR[nrow(Dist)] = Dmin
     
   ntmp = n[imin]+n[jmin]
   mtmp = (n[imin]*m[,imin]+n[jmin]*m[,jmin])/ntmp
   vtmp = (n[imin]*v[,imin] + n[jmin]*v[,jmin] + n[imin]*(m[,imin]-mtmp)^2 + n[jmin]*(m[,jmin]-mtmp)^2 ) / ntmp
      
   Dist <- Dist[-c(imin,jmin), -c(imin,jmin)]
   
   m = m[,-c(imin,jmin)]
   v = v[,-c(imin,jmin)]
   n = n[-c(imin,jmin)]
      
   # update Dist N Nplus LogL
   m = cbind(m , mtmp)  
   v = cbind(v, vtmp)
   n = c(n, ntmp)  
   
   Dtmp = rep(Inf,times=nrow(as.matrix(Dist)))
   
   for (k  in (1:nrow(as.matrix(Dist)))){
   
      ybar      = (n[k]*m[,k]+ntmp*mtmp)/(n[k]+ntmp)
      varpool   = sum(n[k]*v[,k] + ntmp*vtmp + n[k]*(m[,k]-ybar)^2 + ntmp*(mtmp-ybar)^2  ) / (n[k]+ntmp)
      Dtmp[k]   = sum(-n[k]*log(v[,k])-ntmp*log(vtmp))+ (n[k]+ntmp)*log(varpool)
      }
   
   
      Dist = rbind(cbind(Dist,Dtmp),rep(Inf,times=(nrow(as.matrix(Dist))+1)))
      
}
}
phi0= list(mu=m,sigma=sqrt(v),prop=n/sum(n))

invisible(phi0)
}
