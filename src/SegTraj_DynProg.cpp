#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' DynProg_algo_cpp
//'
//' This function finds the best segmentation given a Cost Matrix using a
//' dynamic programming algorithm. C++ implementation of \link{DynProg}
//'
//' @param matD Cost Matrix
//' @param Kmax number of segments
//' @export
// [[Rcpp::export]]

List DynProg_algo_cpp(arma::mat matD, int Kmax) {
  int n = matD.n_rows;
  // I<-matrix(Inf,Kmax,N)
  arma::mat I(Kmax,n);
  std::fill(I.begin(), I.end(), arma::datum::inf);

  // t<-matrix(0,Kmax,N)
  arma::mat tmat(Kmax,n);
  std::fill(tmat.begin(), tmat.end(), 0);

  // I[1,]=matD[1,]
  I.row(0) = matD.row(0);

  // matD=t(matD) ## a quoi ca sert ?
  matD = matD.t();

  // Rcpp::Rcout << "1" << std::endl;
  if(Kmax > 2){
    for(int k = 2; k <= Kmax-1; ++k){
      // Rcpp::Rcout << "k" << k << std::endl;
      for(int L = k; L <= n; ++L){
        mat tmp = I(span(k-2,k-2),span(0,L-2)) + matD(span(L-1,L-1),span(1,L-1));
        uword index_min = tmp.index_min();
        // Rcpp::Rcout << "L" << L << std::endl;

        vec value_min(1);
        std::fill(value_min.begin(),value_min.end(),tmp(index_min));
        I(k-1,L-1) = value_min(0);
        if( ! value_min.is_finite()){
          tmat(k-2,L-1) = arma::datum::inf;
        } else {
          tmat(k-2,L-1) = index_min+1;
        }
        // return(value_min);
        // I[k,L]<-min(I[(k-1),1:(L-1)]+matD[L,2:L])
        //
          //       if(I[k,L]!=Inf){
            //         t[(k-1),L]<-which.min(I[(k-1),1:(L-1)]+matD[L,2:L])
            //       } else {
              //         t[(k-1),L] = Inf
              //       } # end else
      }
    }
  }
  //   I[Kmax,N]<-min(I[(Kmax-1),1:(N-1)]+matD[N,2:N])
  mat tmp2 = I(span(Kmax-2,Kmax-2),span(0,n-2)) + matD(span(n-1,n-1),span(1,n-1));
  uword index_min2 = tmp2.index_min();
  vec value_min2(1);
  std::fill(value_min2.begin(),value_min2.end(),tmp2(index_min2));
  I(Kmax-1,n-1) = value_min2(0);
  //   if(I[Kmax,N]!=Inf){
    //     t[(Kmax-1),N]<-which.min(I[(Kmax-1),1:(N-1)]+matD[N,2:N])
    //   } #end if
  // Rcpp::Rcout << "value min2 = " <<   value_min2  << std::endl;

  // Rcpp::Rcout << "value min2 is finite = " <<   value_min2.is_finite() << std::endl;
  // return I;
  if( value_min2.is_finite() ){
    tmat(Kmax-2,n-1) = index_min2+1;
    // Rcpp::Rcout << "index_min2 = " << index_min2 << std::endl;
  }

  //   t.est<-matrix(0,Kmax,Kmax)
  mat t_est(Kmax,Kmax);
  std::fill(t_est.begin(), t_est.end(), 0);
  vec tmpvec(Kmax);  tmpvec.fill(n);
  //   diag(t.est)<-N
  t_est.diag() = tmpvec;

  // return(tmat);
  for(int K = 2; K <= Kmax; ++K){
    // Rcpp::Rcout << "K= " << K << std::endl;

    for(int k = K-1; k >= 1; k--){
      // Rcpp::Rcout << "k= " << K << std::endl;
      if(t_est(K-1,k) != 0){
        t_est(K-1,k-1) = tmat(k-1,t_est(K-1,k)-1);
      }
    }
  }

  vec J_est = I.col(n-1);
  List ret;
  ret["J.est"] = wrap(J_est);
  ret["t.est"] = wrap(t_est);
  return ret;
  //   for (K in 2:Kmax){
    //     for (k in seq(K-1,1,by=-1)){
      //       if(t.est[K,k+1]!=0){
        //         t.est[K,k]<-t[k,t.est[K,k+1]]
        //       } #end if
      //     } #end k
    //   } #end K
  //     t.est[which(is.na(t.est))]=Inf
  //     list(J.est = I[,N],t.est = t.est)

  //
    // #cat("Kmax", Kmax,"\n")
}
