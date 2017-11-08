#' Test function generating fake data
#' @export
#'
test_data <- function() {
  xcl1 = stats::rnorm(90, mean = 0, sd = 1)
  xcl2 = stats::rnorm(90, mean = 10, sd = 1)
  xcl3 = stats::rnorm(90, mean = 20, sd = 1)
  ycl1 = stats::rnorm(90, mean = 0, sd = 1)
  ycl2 = stats::rnorm(90, mean = 10, sd = 1)
  ycl3 = stats::rnorm(90, mean = 20, sd = 1)
  data = data.frame(x = c(xcl1[1:30],
                          xcl2[1:30],
                          xcl3[1:30],
                          xcl1[31:60],
                          xcl2[31:60],
                          xcl3[31:60],
                          xcl1[61:90],
                          xcl2[61:90],
                          xcl3[61:90]),
                    y = c(ycl1[1:30],
                          ycl2[1:30],
                          ycl3[1:30],
                          ycl1[31:60],
                          ycl2[31:60],
                          ycl3[31:60],
                          ycl1[61:90],
                          ycl2[61:90],
                          ycl3[61:90]))
  true_data_segclust = rep(rep(seq(1,3),each = 30),3)
  true_data_seg = rep(seq(1,9),each = 30)
  return( list( data = data, true_data_seg = true_data_seg, true_data_segclust = true_data_segclust))
}
