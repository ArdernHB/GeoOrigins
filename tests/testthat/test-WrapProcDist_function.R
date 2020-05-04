context("test-WrapProcDist_function")

test_that("Distance table populates correctly",{
  exampleRun <- ProcDistanceTable(Rpraetor$LMs)


  expect_equal(dim(exampleRun)[1], dim(exampleRun)[2])
  expect_equal(dim(exampleRun)[1], dim(Rpraetor$LMs)[3])
  expect_equal(sum(diag(exampleRun)), 0)
  expect_equal(sum(exampleRun[,1]), sum(exampleRun[1,]))


  skip('Parallel processing function testing been skipped as the use of multiple cores can cause testing problems')
  exampleRunPar <- ProcDistanceTablePar(Rpraetor$LMs)
  expect_equal(sum(exampleRun[,1]), sum(exampleRunPar[,1]))
  expect_equal(sum(diag(exampleRun)), sum(diag(exampleRunPar)))
  expect_equal(dim(exampleRunPar)[1], dim(Rpraetor$LMs)[3])
})




