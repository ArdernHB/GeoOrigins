context("test-Geographic_distance_functions")

test_that("multiplication works", {

  ExampleRun <- GeoDist2Point(RefLatLongs = Rpraetor$Lat.Long[-1,],TargetLatLong = Rpraetor$Lat.Long[1,])
  ExampleRunTable <- GeoDist2PointTable(RefLatLongs = Rpraetor$Lat.Long, IDs = rownames(Rpraetor$Lat.Long))

  expect_equal(ExampleRun, ExampleRunTable[1,-1])
  expect_equal(dim(ExampleRunTable)[1], dim(ExampleRunTable)[2])
  expect_equal(dim(ExampleRunTable)[1]-1, length(ExampleRun))
  expect_equal(rownames(ExampleRunTable), colnames(ExampleRunTable))

})
