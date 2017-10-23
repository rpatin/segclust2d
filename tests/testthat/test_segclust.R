context("Segmentation and segclust Testing")

test_that("Test of segmentation", {
  datalist = test_data()
  data <- datalist[["data"]]
  true_data_seg <- datalist[["true_data_seg"]]
  test_seg <- segmentation(data,Kmax = 10,lmin=5,scale.variable = T,type = "home-range",seg.var = c("x","y"))
  test_res <- augment(test_seg)$state
  expect_equal(test_res, true_data_seg)
})

test_that("Test of segmentation/clustering", {
  datalist = test_data()
  data <- datalist[["data"]]
  true_data_segclust <- datalist[["true_data_segclust"]]
  test_segclust <- segclust(data,Kmax = 10,lmin=5,ncluster = 2:5,scale.variable = T,type = "behavior",seg.var = c("x","y"))
  test_res <- augment(test_segclust)$state
  expect_equal(test_res, true_data_segclust)
})
