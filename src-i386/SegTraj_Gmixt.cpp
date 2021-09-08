// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//' arma_repmat
//'
//' C++ Armadillo version for repmat function. Repeat a matrix in bloc.
//'
//' @param A matrix
//' @param n number of repetition in line
//' @param m number of repetition in column
// [[Rcpp::export]]

arma::mat arma_repmat(arma::mat A, int n, int m){
  return arma::repmat(A,n,m);
}

//' Gmixt_algo_cpp
//'
//' Internal C++ algorithm for computing the cost matrix.
//'
//' @param zi vector of observations
//' @param lgi vector of indices
//' @param P number of class
//' @param mvec vector of means for each class
//' @param wk temporary vector for calculations
//' @param svec vector of standard deviations for each class
//' @param prop mixture vector
// [[Rcpp::export]]

arma::mat Gmixt_algo_cpp(arma::vec zi, arma::vec lgi,
                   int P, arma::vec mvec, arma::vec wk,
                   arma::vec svec, arma::vec prop) {
  // zi.print();
  // lgi.print();
  arma::mat tmp = zi/lgi;
  arma::mat tmp_repmat = repmat(tmp.t(),P,1);
  tmp_repmat.each_col() -= mvec;
  arma::mat dkp = square(tmp_repmat);
  dkp.each_row() += wk.t();
  dkp.each_col() /= square(svec);
  dkp.each_col() += log(2*datum::pi*square(svec));
  dkp.each_row() %= (-0.5 * lgi.t());
  dkp.each_col() += log(prop);
  arma::mat Amax = max(dkp,0);
  arma::mat Aprov = dkp;
  Aprov.each_row() -= Amax;
  Aprov = exp(Aprov);
  arma::mat Asum = sum(Aprov,0);
  return -log(Asum) - Amax;
}


//' Gmixt_simultanee_fullcpp
//'
//' C++ function replacing \link{Gmixt_simultanee}
//'
//' @param Don Bivariate Signal
//' @param lmin minimum length of segments
//' @param prop mixture parameters
//' @param mu mean parameters
//' @param s standard deviation parameters
// [[Rcpp::export]]

arma::mat Gmixt_simultanee_fullcpp(arma::mat Don,int lmin, arma::rowvec prop, arma::mat mu, arma::mat s){

  //   prop = phi$prop
  // P = length(phi$prop)
  int P = prop.size();
  //   m    = phi$mu
  //   s    = phi$sigma
  //   n = dim(Don)[2]
  int n = Don.n_cols;
  //   G = list()
  arma::cube G(n,n,2);
  std::fill(G.begin(), G.end(), datum::inf);
  //
  // #rappel: on fait du plus court chemin
  // #        donc on prend -LV
  // # G[[signal]][1, l] contains the likelihood of the segment starting in 1 and finishing
  //
  //     G[[signal]] <- arma::matrix(Inf,ncol=n,nrow=n)
  //   for (signal in 1:2){
  IntegerVector seqn2 = seq_len(n);
  arma::vec seqn = as<arma::vec>(seqn2);
  // Rcpp::Rcout << "1" << std::endl;
  // seqn.print(" seqn ");
  seqn.shed_rows(0,lmin-2);
  // Rcpp::Rcout << "2" << std::endl;

  for(int signal = 0; signal <= 1; signal++){
    //     z = Don[signal,]
    // Rcpp::Rcout << "signal = " << signal << std::endl;
    arma::vec z = Don.row(signal).t();
    //     lg  = lmin:n  # possible position for the end of  first segment
    //       zi  = cumsum(z) # z cumul
    arma::vec zi  = cumsum(z);
    //       zi=zi[lg]
    zi.shed_rows(0,lmin-2);
    // arma::vec zi_sub = zi.subvec(lmin-1,n-1);
    //     z2  = z^2
    arma::vec z2 = square(z);
    //     z2i = cumsum(z2)
    arma::vec z2i = cumsum(z2);
    //       z2i=z2i[lg]
    // arma::vec z2i_sub = z2i.subvec(lmin-1,n-1);
    z2i.shed_rows(0,lmin-2);
    // return seqn;

    //       wk  <- (z2i/lg-(zi/lg)^2)
    arma::vec wk = z2i/seqn - square(zi/seqn);

    // Rcpp::Rcout << "1" << std::endl;
    // G[[signal]][1,lmin:n] = -log(apply(Aprov,2,sum)) - A_max
    arma::mat tmp = Gmixt_algo_cpp(zi, seqn, P, mu.row(signal).t(), wk, s.row(signal).t(), prop.t());
    arma::cube tmp2((const double*)tmp.begin(), 1, n-lmin+1,1);
    G(span(0,0),span(lmin-1,n-1),span(signal,signal)) = tmp2;
    // return tmp2;
    //         for (i in (2:(n-lmin+1))) {
    for (int i = 2; i <= (n-lmin+1); ++i)  {
      // Rcout << "loop " << i << std::endl;

      //           ni  = n-i-lmin+3
      //           z2i = z2i[2:ni]-z2[i-1]

      z2i.shed_row(0);
      z2i -= z2(i-2);
      // Rcout << "z2i =  " << z2i.size() << std::endl;

      //           zi  = zi[2:ni]-z[i-1]
      zi.shed_row(0);
      zi -= z(i-2);
      //           lgi = lmin:(n-i+1)
      arma::vec lgi = seqn.subvec(0,n-i-lmin+1);

      //           wk = (z2i)/(lgi)-(zi/(lgi))^2
      arma::vec wki = z2i/lgi - square(zi/lgi);

      arma::mat tmpi = Gmixt_algo_cpp(zi, lgi, P, mu.row(signal).t(), wki, s.row(signal).t(), prop.t());

      arma::cube tmp2((const double*)tmpi.begin(), 1, n-i-lmin+2,1);
      // return tmp2;

      G(span(i-1,i-1),span(i+lmin-2,n-1),span(signal,signal)) = tmp2;
    }
  }
  arma::mat G2 = sum(G,2);
  return G2;
}

//' cumsum_cpp
//'
//' C++ function for cumulative sum (replacing R cumsum)
//'
//' @param x Numerical Vector
// [[Rcpp::export]]

NumericVector cumsum_cpp(NumericVector x){
  return(cumsum(x));
}