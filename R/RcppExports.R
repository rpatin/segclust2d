# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' DynProg_algo_cpp
#'
#' This function finds the best segmentation given a Cost Matrix using a
#' dynamic programming algorithm. C++ implementation of \link{DynProg}
#'
#' @param matD Cost Matrix
#' @param Kmax number of segments
#' @export
DynProg_algo_cpp <- function(matD, Kmax) {
    .Call(`_segclust2d_DynProg_algo_cpp`, matD, Kmax)
}

#' logdens_simultanee_cpp
#'
#' Calculate logdensity of a bivariate signal
#'
#' @rdname logdens_simultanee
#' @param mu mean parameter for each signal
#' @param sigma standard deviation parameter for each signal
#' @param prop mixture parameter
logdens_simultanee_cpp <- function(xk, mu, sigma, prop) {
    .Call(`_segclust2d_logdens_simultanee_cpp`, xk, mu, sigma, prop)
}

#' apply_rowSums
#'
#' Internal function for Expectation-Maximization (EM) algorithm.
#'
#' @param rupt current estimated breaks in signal
#' @param x bivariate signal
apply_rowSums <- function(rupt, x) {
    .Call(`_segclust2d_apply_rowSums`, rupt, x)
}

#' colsums_sapply
#'
#' Internal function for Expectation-Maximization (EM) algorithm.
#'
#' @param rupt current estimated breaks in signal
#' @param x bivariate signal
#' @param mu mean parameter for each signal
#' @param tau tau
#' @param i number of signal
colsums_sapply <- function(i, rupt, x, mu, tau) {
    .Call(`_segclust2d_colsums_sapply`, i, rupt, x, mu, tau)
}

#' arma_repmat
#'
#' C++ Armadillo version for repmat function. Repeat a matrix in bloc.
#'
#' @param A matrix
#' @param n number of repetition in line
#' @param m number of repetition in column
arma_repmat <- function(A, n, m) {
    .Call(`_segclust2d_arma_repmat`, A, n, m)
}

#' Gmixt_algo_cpp
#'
#' Internal C++ algorithm for computing the cost matrix.
#'
#' @param zi vector of observations
#' @param lgi vector of indices
#' @param P number of class
#' @param mvec vector of means for each class
#' @param wk temporary vector for calculations
#' @param svec vector of standard deviations for each class
#' @param prop mixture vector
Gmixt_algo_cpp <- function(zi, lgi, P, mvec, wk, svec, prop) {
    .Call(`_segclust2d_Gmixt_algo_cpp`, zi, lgi, P, mvec, wk, svec, prop)
}

#' Gmixt_simultanee_fullcpp
#'
#' C++ function replacing \link{Gmixt_simultanee}
#'
#' @param Don Bivariate Signal
#' @param lmin minimum length of segments
#' @param prop mixture parameters
#' @param mu mean parameters
#' @param s standard deviation parameters
Gmixt_simultanee_fullcpp <- function(Don, lmin, prop, mu, s) {
    .Call(`_segclust2d_Gmixt_simultanee_fullcpp`, Don, lmin, prop, mu, s)
}

#' cumsum_cpp
#'
#' C++ function for cumulative sum (replacing R cumsum)
#'
#' @param x Numerical Vector
cumsum_cpp <- function(x) {
    .Call(`_segclust2d_cumsum_cpp`, x)
}

