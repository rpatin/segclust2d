# Gmean_simultanee
#' Gmean_simultanee  calculates the cost matrix for a segmentation model with changes in the mean and same variance for all signals
#' @param Don the bivariate signal
#' @param lmin minimum size for a segment, default value is 2
#' @param sameVar TRUE if segments have the same variance
#' @return the cost matrix G(i,j) which contains the variance of the data between point (i+1) to point j
#' @export
Gmean_simultanee<-function(Don,lmin,sameVar=FALSE)
{

  n = dim(Don)[2]

  if(sameVar){
    ## every element of the list is the cost motrix for one signal
    matD_list=lapply(1:2,function(nb) {
      Res=matrix(Inf,n,n)
      z=Don[nb,] #depend de la forme des donnees
      z2=z^2
      z2i=cumsum(z2)
      zi=cumsum(z)
      z2i=z2i[lmin:n]
      zi=zi[lmin:n]
      Res[1,lmin:n]=z2i-((zi^2)/(lmin:n))
      nl=n-lmin+1
      for (i in 2:nl)
      {
        ni=n-i-lmin+3
        z2i=z2i[2:ni]-z2[i-1]
        zi=zi[2:ni]-z[i-1]
        deno<-((i+lmin-1):n)-i+1
        Res[i,(i+lmin-1):n]=z2i-((zi^2)/deno)
      }
      return(Res)
    })
  } else {
    if (lmin<5){
      cat("lmin =",lmin," and should be > 5 when sameV =FALSE to avoid variance estimation instability  ", "\n")
    }
    # segmentation with heterogeneous variances

    matD_list=lapply(1:2,function(nb) {
      Res=matrix(Inf,n,n)
      z=Don[nb,] #depend de la forme des donnees
      z2=z^2
      z2i=cumsum(z2)
      zi=cumsum(z)
      z2i=z2i[lmin:n]
      zi=zi[lmin:n]

      Res[1,lmin:n]=(lmin:n)*log(( z2i-(zi^2)/(lmin:n) ) / (lmin:n) )
      nl=n-lmin+1
      for (i in 2:nl){
        ni=n-i-lmin+3
        z2i=z2i[2:ni]-z2[i-1]
        zi=zi[2:ni]-z[i-1]
        Res[i,(i+lmin-1):n]=(lmin:(n-i+1))*(log((z2i-(zi^2)/(lmin:(n-i+1)))/(lmin:(n-i+1))))

      }
      return(Res)
    })
  }
  matD=Reduce("+",matD_list)
  invisible(matD)
}
