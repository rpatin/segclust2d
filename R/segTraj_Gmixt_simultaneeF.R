# calculates the cost matrix for a segmentation/clustering model
# lmin : minimum size for a segment
# phi  : parameters of the mixture
# P    : number of clusters
# G(i,j) = mixture density for segment between points (i+1) to j
#        = \sum_{p=1}^P \log (\pi_p f(y^{ij};\theta_p))
# Rq: this density if factorized in order to avoid numerical zeros in the log
Gmixt_simultaneeF <- function(Don,lmin,phi,P){

m    = phi$mu
s    = phi$sigma
prop = phi$prop

n = dim(Don)[2]
#x = x'


GG = matrix(0,ncol=n,nrow=n)


for (signal in 1:2){
x=Don[signal,]
xi  = cumsum(x) 
lg  = lmin:n
xi  = xi[lg]
x2  = x^2
x2i = cumsum(x2)
x2i = x2i[lg]

wkF  = repmat(t( x2i/lg-(xi/lg)^2 ),P,1)

#wk=repmat(wk,P,1)

dkpF   = (repmat(t(xi),P,1)/repmat(t(lg),P,1)-repmat(m[signal,],1,n-1))^2
AF     = (wkF+dkpF)/repmat(s[signal,]^2,1,n-lmin+1)+log(2*pi*repmat(s[signal,]^2,1,n-1))
AF     = -0.5*repmat(t(lg),P,1)*AF +(repmat(log(prop),1,n-1))
AF_max = apply(AF,2,max)
AF     = exp(AF-repmat(t (AF_max) ,P,1))


#rappel: on fait du plus court chemin
#        donc on prend -LV
GG[1,lmin:n] = -log(apply(AF,2,sum)) - AF_max

for (i in (2:(n-lmin+1))) {
   ni  = n-i-lmin+3
   x2i = x2i[2:ni]-x2[i-1]
   xi  = xi[2:ni]-x[i-1]
   lgi = lmin:(n-i+1)
   wkF  = repmat(t(x2i)/(lgi)-(xi/(lgi))^2,P,1)
   dkpF = (repmat(t(xi),P,1)/repmat(t(lgi),P,1)-repmat(m[signal,],1,ni-1))^2
   AF   = (wkF+dkpF)/repmat(s[signal,]^2,1,ni-1)+log(2*pi*repmat(s[signal,]^2,1,ni-1))
   AF   = -0.5*repmat(t(lgi),P,1)*AF +(repmat(log(prop),1,ni-1))
   AF_max = apply(AF,2,max)
   AF     = exp(AF-repmat(t (AF_max) ,P,1))

   GG[i,(i+lmin-1):n] =  GG[i,(i+lmin-1):n]-log(apply(AF,2,sum)) - AF_max
}}
for (i in (lmin-1):n){
  for (j in 1:i){
  GG[i,j]=Inf
  }}
invisible(GG)
}
