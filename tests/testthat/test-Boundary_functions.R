context("test-Boundary_functions")

test_that("BoundaryFinder produces correct outputs",{
  testSelect <- c(1,5,20,30,35)
  RatDistMat <- ProcDistanceTable(Rpraetor$LMs[,,testSelect])

  Long.Range <- c(floor(min(Rpraetor$Lat.Long$Long[testSelect])),
                  ceiling(max(Rpraetor$Lat.Long$Long[testSelect])))
  Lat.Range <- c(floor(min(Rpraetor$Lat.Long$Lat[testSelect])),
                 ceiling(max(Rpraetor$Lat.Long$Lat[testSelect])))

  rThres <- IDbyDistanceDistInputCCV(LatLongs = Rpraetor$Lat.Long[testSelect,],
                                     DistDataMat = RatDistMat,
                                     Verbose = TRUE,
                                     ProvConfidence = .95,
                                     PrintProg = FALSE,
                                     Method = 'Spearman')

  LowResRsamp.Rp <- c(5,5)
  Boundaryfinding <- BoundaryFinder(LatLongs = Rpraetor$Lat.Long[testSelect,],
                                    RefDistMat = RatDistMat,
                                    LongRange = Long.Range,
                                    LatRange = Lat.Range,
                                    RangeSamp = LowResRsamp.Rp,
                                    PlotValCor = rThres$`Provenancing.Correlation.95%.Confidence`,
                                    ExpandMap = c(0,0),
                                    RefIDs = rownames(Rpraetor$Lat.Long[testSelect,]),
                                    DataDump = FALSE,
                                    IgnorePrompts = TRUE,
                                    Method = 'Spearman')



  expect_output(str(Boundaryfinding), "List of 2")
  expect_equal(dim(Boundaryfinding[[1]])[3], dim(Rpraetor$LMs[,,testSelect])[3])

  testSelect2 <- 20:35

  SpecimenLoc <- cbind(chr2nu(Rpraetor$Lat.Long$Long[testSelect2]), chr2nu(Rpraetor$Lat.Long$Lat[testSelect2]))
  DistPol <-  sp::Polygon(SpecimenLoc[grDevices::chull(SpecimenLoc),])

  expect_s4_class(DistPol, "Polygon")

  DistPols <-  sp::Polygons(list(DistPol),1)

  expect_s4_class(DistPols, "Polygons")

  DistSpatPol <-  sp::SpatialPolygons(list(DistPols))

  expect_s4_class(DistSpatPol, 'SpatialPolygons')


  expect_output(sp::proj4string(DistSpatPol), NA)

  expect_warning(raster::area(DistSpatPol))

  sp::proj4string(DistSpatPol) <- "+proj=longlat +ellps=WGS84"
  expect_equal(raster::area(DistSpatPol), 114824797235)



})

