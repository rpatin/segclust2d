// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

arma::mat arma_repmat(arma::mat A, int n, int m){
  return arma::repmat(A,n,m);
}


// [[Rcpp::export]]

arma::mat arma_repmat_transpose_divide(arma::vec x, arma::vec y, int n, int m){
  arma::mat tmp = x/y;
  return arma::repmat(tmp.t(),n,m);
}



// [[Rcpp::export]]
arma::mat sweep_row_plus(arma::mat dkp, arma::vec wk) {
  return dkp.each_row() +  wk.t();
}

// [[Rcpp::export]]
arma::mat sweep_col_divide(arma::mat A, arma::vec s) {
  return A.each_col() / s;
}

// [[Rcpp::export]]
arma::mat sweep_col_plus(arma::mat A, arma::vec s) {
  return A.each_col() + s;
}

// [[Rcpp::export]]
arma::mat sweep_row_times(arma::mat A, arma::vec s) {
  return A.each_row() % s.t();
}

// [[Rcpp::export]]
arma::mat sweep_row_times2(arma::mat A, arma::vec s) {
  int A_nrow = A.n_rows;
  for(int i = 0; i < A_nrow; i++){
    A.row(i) = A.row(i)*s*(-0.5);
  }
  return A;
}

// [[Rcpp::export]]
arma::rowvec apply_col_max(arma::mat A) {
  return max(A, 0);
}

// [[Rcpp::export]]
arma::rowvec apply_col_sum(arma::mat A) {
  return sum(A, 0);
}

// [[Rcpp::export]]
arma::mat Gmixt_algo_cpp(arma::vec zi, arma::vec lgi,
                         int P, arma::vec mvec, arma::vec wk,
                         arma::vec svec, arma::vec prop) {
  arma::mat tmp = zi/lgi;
  arma::mat tmp_repmat = arma::repmat(tmp.t(),P,1);
  tmp_repmat.each_col() -= mvec;
  arma::mat dkp = square(tmp_repmat);
  dkp.each_row() += wk.t();
  dkp.each_col() /= square(svec);
  dkp.each_col() += log(2*arma::datum::pi*square(svec));
  dkp.each_row() %= (-0.5 * lgi.t());
  dkp.each_col() += log(prop);
  arma::mat Amax = max(dkp,0);
  arma::mat Aprov = dkp;
  Aprov.each_row() -= Amax;
  Aprov = exp(Aprov);
  arma::mat Asum = sum(Aprov,0);
  return -log(Asum) - Amax;
}
