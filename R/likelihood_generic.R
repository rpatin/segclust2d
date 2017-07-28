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

#' Generic function for mBIC
#' @export

mBIC <- function (x, ...) {
  UseMethod("BIC", x)
}

#' Default function for mBIC
#' @export

mBIC.default <- function (x, ...) {
  message("No default function for mBIC")
}

#' Generic function for lavielle_penalty
#' @export

lavielle_penalty <- function (x, ...) {
  UseMethod("lavielle_penalty", x)
}

#' Default function for mBIC
#' @export

lavielle_penalty.default <- function (x, ...) {
  message("No default function for lavielle_penalty")
}
