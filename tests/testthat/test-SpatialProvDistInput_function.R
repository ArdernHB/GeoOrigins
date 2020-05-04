context("test-SpatialProvDistInput_function")

test_that("Correlation output is as expected", {

  ExampleDist <- ProcDistanceTable(Rpraetor$LMs)
  GeographicDist <- GeoDist2Point(RefLatLongs = Rpraetor$Lat.Long[-1,], TargetLatLong = Rpraetor$Lat.Long[1,])

  CorResSpearman <- suppressWarnings(stats::cor.test(x = ExampleDist[1,-1], y = GeographicDist, method = 'spearman'))
  expect_s3_class(CorResSpearman, 'htest')
  expect_equal(length(CorResSpearman),8)

  CorResPearson <- suppressWarnings(stats::cor.test(x = ExampleDist[1,-1], y = GeographicDist, method = 'pearson'))
  expect_s3_class(CorResPearson, 'htest')
  expect_equal(length(CorResPearson),9)

})

