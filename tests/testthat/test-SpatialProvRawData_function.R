context("test-SpatialProvRawData_function")

test_that("Outputs from external functions are as expected", {


  gpaRes <- shapes::procGPA(Rpraetor$LMs)
  expect_equal(length(gpaRes), 17)
  expect_equal(dim(gpaRes$rotated), dim(Rpraetor$LMs))

  ExampleDist <- shapes::procdist(gpaRes$rotated[,,1], gpaRes$rotated[,,2])
  expect_equal(is.numeric(ExampleDist), TRUE)


})
