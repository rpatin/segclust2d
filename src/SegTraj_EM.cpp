#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' logdens_simultanee_cpp
//'
//' Calculate logdensity of a bivariate signal
//'
//' @rdname logdens_simultanee
//' @param mu mean parameter for each signal
//' @param sigma standard deviation parameter for each signal
//' @param prop mixture parameter
// [[Rcpp::export]]

NumericVector logdens_simultanee_cpp(arma::mat xk, arma::mat mu, arma::mat sigma, arma::vec prop) {
  int P = prop.size();
  int nk = xk.n_cols;
  arma::vec tmp(P);


  for( int i = 0; i < P; ++i){

    // Rcpp::Rcout << "sigma(1,i)" << sigma(1,i) << std::endl;
    // Rcpp::Rcout << "xk.row(1)" << xk.row(1)<< std::endl;
    // Rcpp::Rcout << "mu(1,i)" << mu(1,i)<< std::endl;
    // Rcpp::Rcout << "1" << << std::endl;
    // Rcpp::Rcout << "1" << << std::endl;

    tmp(i) = - nk * log( sqrt(2*datum::pi)*sigma(0,i) ) -
      0.5 * sum( square(xk.row(0) - mu(0,i)) ) / pow(sigma(0,i),2.0) -
      nk * log( sqrt(2*datum::pi)*sigma(1,i) ) -
      0.5 * sum( square(xk.row(1) - mu(1,i)) ) / pow(sigma(1,i),2.0);
    //     -nk*log(sqrt(2*pi)*s[1,p])
    // - 0.5*sum ( (xk[1,]  - m[1,p])^2  )/s[1,p]^2
    // -nk*log(sqrt(2*pi)*s[2,p])
    // - 0.5*sum (    (xk[2,]  - m[2,p])^2  )/s[2,p]^2

  }
  return wrap(tmp);
}

//' apply_rowSums
//'
//' Internal function for Expectation-Maximization (EM) algorithm.
//'
//' @param rupt current estimated breaks in signal
//' @param x bivariate signal
// [[Rcpp::export]]

arma::mat apply_rowSums(arma::mat rupt, arma::mat x){
  int nrupt = rupt.n_rows;
  arma::mat final(2,nrupt);
  for(int r = 0; r < nrupt; ++r){
    // arma::mat tmp =  sum(x.cols(rupt(r,0), rupt(r,1));
    // Rcpp::Rcout << "tmp" << tmp.n_rows << std::endl;
    // Rcpp::Rcout << "tmp" << tmp.n_cols << std::endl;
    final.col(r) = sum(x.cols(rupt(r,0)-1, rupt(r,1)-1),1);
  }
  return final;
}


//' colsums_sapply
//'
//' Internal function for Expectation-Maximization (EM) algorithm.
//'
//' @param rupt current estimated breaks in signal
//' @param x bivariate signal
//' @param mu mean parameter for each signal
//' @param tau tau
//' @param i number of signal
// [[Rcpp::export]]


arma::mat colsums_sapply(int i, arma::mat rupt, arma::mat x, arma::mat mu, arma::mat tau){
  // s[i,] =  colSums( tau*(
  //   sapply(1:P, function(p) {
  //     apply(rupt,1,FUN=function(y) sum((x[i,y[1]:y[2]]-m[i,p])^2   ))
  //   })#sapply
  // )#tau
  // ) #colsums
  int P = tau.n_cols;
  int nrupt = rupt.n_rows;
  arma::mat sapply_mat(nrupt,P);
  // Rcpp::Rcout << "  x(span(i-1,i-1),span(rupt(r,0)-1, rupt(r,1)-1))" <<   x(span(i-1,i-1),span(rupt(0,0)-1, rupt(0,1)-1)) << std::endl;
  for(int p = 0; p < P; ++p){
    // Rcpp::Rcout << "p = " << p << std::endl;
    for(int r = 0; r < nrupt; ++r){
      // Rcpp::Rcout << "r = " << r << std::endl;
      // Rcpp::Rcout << "mu(i-1,r-1) = " << mu(i-1,p) << std::endl;
      // return x(span(i-1,i-1),span(rupt(r,0)-1, rupt(r,1)-1)) - mu(i-1,p);
      sapply_mat(r,p) =  accu( pow(x(span(i-1,i-1),span(rupt(r,0)-1, rupt(r,1)-1)) - mu(i-1,p),2) );
      // return accu( pow(x(span(i-1,i-1),span(rupt(r,0)-1, rupt(r,1)-1)) - mu(i-1,p),2) );
    }
  }
  // return sapply_mat;
  return sum(sapply_mat % tau,0);
}
