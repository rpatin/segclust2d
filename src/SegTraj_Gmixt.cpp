// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::mat arma_repmat(arma::mat A, int n, int m){
  return arma::repmat(A,n,m);
}


// [[Rcpp::export]]
Rcpp::NumericMatrix sweep_row_plus(Rcpp::NumericMatrix dkp, Rcpp:: NumericVector wk) {
  int dkp_nrow = dkp.nrow();
  int dkp_ncol = dkp.ncol();
  for(int i = 0; i < dkp_nrow; i++){
    for(int j = 0; j < dkp_ncol; j++){
      dkp(i,j) += wk(j);
    }
  }
  return dkp;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sweep_col_divide(Rcpp::NumericMatrix A, Rcpp:: NumericVector s) {
  int A_nrow = A.nrow();
  int A_ncol = A.ncol();
  for(int i = 0; i < A_nrow; i++){
    for(int j = 0; j < A_ncol; j++){
      A(i,j) /= s(i);
    }
  }
  return A;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sweep_col_plus(Rcpp::NumericMatrix A, Rcpp::NumericVector s) {
  int A_nrow = A.nrow();
  int A_ncol = A.ncol();
  for(int i = 0; i < A_nrow; i++){
    for(int j = 0; j < A_ncol; j++){
      A(i,j) += s(i);
    }
  }
  return A;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sweep_row_times(Rcpp::NumericMatrix A, Rcpp::NumericVector s) {
  int A_nrow = A.nrow();
  int A_ncol = A.ncol();
  for(int i = 0; i < A_nrow; i++){
    for(int j = 0; j < A_ncol; j++){
      A(i,j) *= s(j);
    }
  }
  return A;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sweep_row_times2(Rcpp::NumericMatrix A, Rcpp::NumericVector s) {
  int A_nrow = A.nrow();
  int A_ncol = A.ncol();
  for(int i = 0; i < A_nrow; i++){
    for(int j = 0; j < A_ncol; j++){
      A(i,j) *= s(j);
      A(i,j) *= -0.5;
    }
  }
  return A;
}

// [[Rcpp::export]]
Rcpp::NumericVector apply_col_max(Rcpp::NumericMatrix A) {
  int A_nrow = A.nrow();
  int A_ncol = A.ncol();
  Rcpp::NumericVector col_max(A_ncol);
  for(int i = 0; i < A_nrow; i++){
    for(int j = 0; j < A_ncol; j++){
      col_max(j) = Rcpp::max(A( Rcpp::_, j));
    }
  }
  return col_max;
}

// [[Rcpp::export]]
Rcpp::NumericVector apply_col_sum(Rcpp::NumericMatrix A) {
  int A_nrow = A.nrow();
  int A_ncol = A.ncol();
  Rcpp::NumericVector col_sum(A_ncol);
  for(int i = 0; i < A_nrow; i++){
    for(int j = 0; j < A_ncol; j++){
      col_sum(j) = Rcpp::sum(A( Rcpp::_, j));
    }
  }
  return col_sum;
}


