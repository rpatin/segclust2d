# repmat
#' repmat repeats a matrix
#' @param a the base matrix
#' @param n number of repetition in lines
#' @param m number of repetition in columns
#' @return a matrix with n repeats of a in lines et m in columns
repmat <- function(a,n,m) {kronecker(matrix(1,n,m),a)} 
