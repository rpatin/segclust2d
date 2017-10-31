#' Generic function for augment
#'
#' see broom::augment for more informations
#' @param x object to be augmented
#' @export

augment <- function (x, ...) {
  UseMethod("augment", x)
}


