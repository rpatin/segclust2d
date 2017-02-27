#' Generic function for likelihood
#' @export

likelihood <- function (x, ...) {
  UseMethod("likelihood", x)
}

#' Default function for likelihood
#' @export

likelihood.default <- function (x, ...) {
  message("No default function for likelihood")
}
