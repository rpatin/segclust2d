#' Internal function for choosing optimal number of segment
#'
#' Internal function for choosing optimal number of segment using Marc Lavielle's method. From Emilie Lebarbier. Method based on identifying breaks in the slope of the contrast.

#' @param J likelihood for each number of segment
#' @param S threshold for choosing the number of segment. See
#'   adehabitatLT::chooseseg
#' @return  a list with optimal number of segment and full data.frame of the
#'   calculus
#'
#' @examples
#' \dontrun{segmentation(data, diag.var=c("dist","angle"),
#' order.var='dist',type='hmm', hmm.model=mod1.hmm)}
chooseseg_lavielle <- function(J, S=0.75)
{
  Kmax = length(J)
  Kseq = 1:Kmax
  Kmax=length(Kseq)
  Jtild=(J[Kmax]-J)/(J[Kmax]-J[1])*(Kseq[Kmax]-Kseq[1])+1
  D=diff(diff(Jtild))
  if (length(which(D>=S))>0){
    Kh <- max(which(D>=S))+1
  }else {
    Kh <- 1
  }
  Kh=Kseq[Kh]
  return(list("Kopt"= Kh, "lavielle"= D))
}

# chooseseg_BIC <- function(J.est, S=0.75)
# {
#   Kmax = length(J.est)
#   J.est2 <- ((J.est[Kmax]-J.est)/(J.est[Kmax]-J.est[1]))*(Kmax-1)+1
#   D <- c(Inf, sapply(2:(length(J.est2)-1), function(K) J.est2[K-1] - 2*J.est2[K] + J.est2[K+1]))
#   df <- data.frame(K=1:(Kmax-1), D=D, J.est = J.est[-Kmax])
#   cons <- c(as.numeric(J.est[-1]<J.est[-length(J.est)]))
#   if (length(df$K[df$D>S&cons==1])>=1) {
#     Kopt <- max(df$K[df$D>S&cons==1])
#   } else {
#     Kopt <- 1
#   }
#   return(list("Kopt"= Kopt, "lavielle"= df))
# }
