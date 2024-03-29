% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R, R/segTraj_logdens_simultanee.R
\name{logdens_simultanee_cpp}
\alias{logdens_simultanee_cpp}
\alias{logdens_simultanee}
\title{logdens_simultanee_cpp}
\usage{
logdens_simultanee_cpp(xk, mu, sigma, prop)

logdens_simultanee(xk, phi)
}
\arguments{
\item{xk}{the bivariate signal}

\item{mu}{mean parameter for each signal}

\item{sigma}{standard deviation parameter for each signal}

\item{prop}{mixture parameter}

\item{phi}{parameters of the mixture, P components}
}
\value{
the value of the log density
}
\description{
Calculate logdensity of a bivariate signal
}
