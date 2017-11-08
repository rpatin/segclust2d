#' Generic function for likelihood
#' @param x object from which likelihood can be extracted
#' @param ... additional parameters
#' @export

likelihood <- function (x, ...) {
  .Deprecated("logLik")
  UseMethod("likelihood", x)
}

