#' Generic function for augment
#' @export

augment <- function (x, ...) {
  UseMethod("augment", x)
}


