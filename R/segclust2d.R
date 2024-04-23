#' segclust2d: tools for segmentation of animal GPS movement data
#'
#' Provides two methods for segmentation and joint segmentation/clustering of
#' bivariate time-series. Originally intended for ecological segmentation
#' (home-range and behavioural modes) but easily applied on other series, the
#' package also provides tools for analysing outputs from R packages moveHMM
#' and marcher.
#'
#' The segmentation method is a bivariate extension of Lavielle's method
#' available in adehabitatLT (Lavielle 1999; and 2005). This method rely on
#' dynamic programming for efficient segmentation.
#'
#' The segmentation/clustering method alternates steps of dynamic programming
#' with an Expectation-Maximization algorithm. This is an extension of Picard et
#' al (2007) method (formerly available in cghseg package) to the bivariate
#' case. 
#' 
#' The full description of the method is published in Patin et al. (2020).
#'
#' References:
#'
#'Lavielle, M. (1999) Detection of multiple changes in a sequence of dependent
#'variables. \emph{Stochastic Processes and their Applications}, \bold{83}:
#'79--102.
#'
#'Lavielle, M. (2005) Using penalized contrasts for the change-point problem.
#'Report number 5339, Institut national de recherche en informatique et en
#'automatique.
#'
#'Patin, R., Etienne, M. P., Lebarbier, E., Chamaille-Jammes, S., 
#'& Benhamou, S. (2020). Identifying stationary phases in
#' multivariate time series for highlighting behavioural modes 
#' and home range settlements. \emph{Journal of Animal Ecology}, 
#' 89(1), 44-56.
#'
#'
#'Picard, F., Robin, S., Lebarbier, E. and Daudin, J.-J. (2007), A
#'Segmentation/Clustering Model for the Analysis of Array CGH Data.
#'\emph{Biometrics}, 63: 758-766. doi:10.1111/j.1541-0420.2006.00729.x
#'
#' @name segclust2d
#' @useDynLib segclust2d, .registration=TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"