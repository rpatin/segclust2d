#' Generic function for likelihood
#' @export

likelihood <- function (x, ...) {
  .Deprecated("logLik")
  UseMethod("likelihood", x)
}


#' Generic function for BIC
#' @export

BIC <- function (x, ...) {
  UseMethod("BIC", x)
}

#' Default function for BIC
#' @export

BIC.default <- function (x, ...) {
  message("No default function for BIC")
}
