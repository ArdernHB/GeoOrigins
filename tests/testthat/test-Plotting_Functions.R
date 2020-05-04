context("test-Plotting_Functions")

test_that("Contour polygon construction works", {

  CorMatrixRes <- Marvalis$Orkney.Boundary$RawCorData[,,1]
  CoordsHeat <- as.data.frame(cbind(expand.grid(as.numeric(rownames(CorMatrixRes)),as.numeric(colnames(CorMatrixRes)))))
  colnames(CoordsHeat) <- c('Lats', 'Longs')

  KSres <- ks::kde(x=CoordsHeat, H=ks::Hpi(x=CoordsHeat), compute.cont=TRUE)
  expect_s3_class(KSres, 'kde')

  ContourExample <- Construct_contour(CoordsHeat)
  expect_equal(length(ContourExample), 3)

})
