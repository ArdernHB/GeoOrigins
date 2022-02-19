


#' Spatially provenance a specimen from distance data adapted for use in parallel processing
#'
#' This function is primarily for compatibility with the parallel processing \code{foreach} package.
#' It takes the distances from an unknown specimen to all reference specimens
#' and uses these distances to calculate a likely spatial provenance. Note that this
#' procedure can only be applied to one unknown specimen at a time. This distance input function
#' allows the user to input any dissimilarity or distance desired and as appropriate for the data.
#' The function has two applications either: calculating a specimens' provenance, or alternatively
#' it can be used to calculate the minimum correlation coefficient needed to correctly identify a
#' known specimen at its true collection location. The second application of this function can work
#' as a correct cross-validation process if looped, but see IDbyDistanceDistInputCCV function which
#' does this automatically in a leave-one-out process. However, if a corss-validation process that
#' removed more than one specimen from the reference dataset at a time is required then it is advised
#' that this be applied using this function.
#' @param DistDataVecPar is a vector of distances from each reference specimen to the specimen of interest (either a specimen with unknown provenance or another reference specimen that the user is interested in validating the provenance of)
#' @param LatLongsPar a matrix of n rows by 2 columns where n is the number of reference specimens in your dataset and the columns are Latitude and Longitude values in that order. These latitude-longitude coordinates should be of the locations of the reference specimens.
#' @param LongRangePar is a vector of 2 elements defining the maximum and minimum Longitude values that the provenancing method should explore. This will also define the mapping range in the final plotted output.
#' @param LatRangePar is a vector of 2 elements defining the maximum and minimum Latitude values that the provenancing method should explore. This will also define the mapping range in the final plotted output.
#' @param MethodPar determines what kind of correlation coefficient should be used, either "Spearman" or "Pearson". Spearman's ranked correlation coefficient does not assume a linear relationship between geographic and trait distances, whereas Pearson's coefficient does.
#' @inheritParams IDbyDistanceRawData
#' @return If Verbose is TRUE then a dataframe of all values for every sampled grid reference is returned. If Verbose is FALSE then only those grid references with the highest correlation values are returned.
#' @details This method also makes use of the \code{\link[stats]{cor.test}} function from the \code{stats} package. When the \code{PrintProg} is set to TRUE, the \code{\link[svMisc]{progress}} function of the \code{svMisc} package is used.
#' The map plotting of this function makes use of the functions of the \code{maps} package.
#'
#' @keywords internal
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @author Ardern Hulme-Beaman



IDbyDistanceDistInputPar <- function(LatLongsPar, DistDataVecPar, LongRangePar, LatRangePar, RangeSamp=10, MethodPar=c('Spearman', 'Pearson')){

  UserInputAssessment(LatLongs=LatLongsPar, DistVec = DistDataVecPar, Method=MethodPar, RefDistMat = NULL, RefData = NULL)

  #making LatLongs a dataframe
  LatLongsPar <- as.data.frame(LatLongsPar)
  colnames(LatLongsPar) <- c("Lats", "Longs")

  ShapeDist <- DistDataVecPar

  #creating an empty object to be populated by results
  CoordsHeat <- NULL


  #this function can come up with correlation values across the entire map
  #here we have the first option which is useful for the validation process
  #by doing this first option for all the specimens (in a loop) we can build a distribution
  #of the correlation values that will correctly cover the specimens true location
  #but see IDbyDistanceRawDataCCV function for looping to be done automatically

  #creating a range that will cover the whole geographic area of interest
  #this is for the function to loop through later on

  if (length(RangeSamp)==1){
    LongSamp <- RangeSamp
    LatSamp <- RangeSamp
  } else if (length(RangeSamp)==2){
    LongSamp <- RangeSamp[2]
    LatSamp <- RangeSamp[1]
  } else if (length(RangeSamp)>2){
    stop("too many dimensions in RangeSamp")
  }

  LongRangeSteps <- (LongRangePar[2]-LongRangePar[1])/(LongSamp-2)
  LatRangeSteps <- (LatRangePar[2]-LatRangePar[1])/(LatSamp-2)

  #this output for Lat/Long ways provides what the loop should sequence through
  Longways <- c(LongRangePar[1]-LongRangeSteps, seq(LongRangePar[1], LongRangePar[2], by = LongRangeSteps), LongRangePar[2]+LongRangeSteps)
  Latways <- c(LatRangePar[1]-LatRangeSteps, seq(LatRangePar[1], LatRangePar[2], by = LatRangeSteps), LatRangePar[2]+LatRangeSteps)

  CorMatrixRes <- matrix(NA, nrow = length(Latways), ncol = length(Longways))
  rownames(CorMatrixRes) <- Latways
  colnames(CorMatrixRes) <- Longways

  GeoDist2PointPar <- function(RefLatLongs, TargetLatLong){
    GeographicDist <- apply(X = RefLatLongs[,2:1], MARGIN = 1, FUN = geosphere::distHaversine, p2 = TargetLatLong[2:1], r=6378.137)
    return(GeographicDist)
  }
  #carrying out iterative analyses across the geographic region defined by LongRange and LatRange
  for (i in 1:length(Longways)){
    for (j in 1:length(Latways)){
      #i <- 1
      #j <- 1

      coord <- c(Latways[j], Longways[i])

      GeographicDist <- GeoDist2PointPar(RefLatLongs = LatLongsPar, TargetLatLong = coord)

      if (MethodPar=='Spearman'){
        #running the correlation to generate r
        CorRes <- suppressWarnings(stats::cor.test(x = ShapeDist, y = GeographicDist, method = "spearman"))
      } else if (MethodPar=='Pearson'){
        CorRes <- suppressWarnings(stats::cor.test(x = ShapeDist, y = GeographicDist, method = "pearson"))
      }

      #populating results matrix
      CorMatrixRes[j,i] <- CorRes$estimate

    }


  }

  CoordsHeat <- CorMatrixRes
  return(CoordsHeat)
}





#' Spatial Provenancing Correct Cross-Validation calculation from distance data using parallel processing
#'
#' This function takes pairwise distances among all reference specimens with known spatial origins
#' and uses those distances to calculate the correllation value that would be required to
#' correctly provenancing them in a test. This is achieved by a leave-one-out procedure. If a
#' cross-validation procedure that removes more than one specimen from the reference dataset is
#' desired then it is recommended that the validate function of the IDbyDistanceDistInput method be used
#' with the validate argument in a loop.
#' @param DistDataMat is a square matrix of pairwise distances among all reference specimens
#' @inheritParams IDbyDistanceRawDataCCV
#' @return If Verbose is set to FALSE then a list of a single object containing the correlation value at the required confidence interval is returned. If Verbose is set to TRUE then a list is returned with two objects: the first is the correlation value at the required confidence interval; the second a dataframe of coordinates and the spatial-trait correlation values at the true locations of each specimen.
#' @details This method also makes use of the \code{\link[stats]{cor.test}} function from the \code{stats} package. When the \code{PrintProg} is set to TRUE, the \code{\link[svMisc]{progress}} function of the \code{svMisc} package is used.
#' The map plotting of this function makes use of the functions of the \code{maps} package.
#'
#' @section Citations:
#'
#' Original S code by Richard A. Becker, Allan R. Wilks. R version by Ray Brownrigg.
#' Enhancements by Thomas P Minka and Alex Deckmyn. (2017). maps: Draw Geographical Maps. R
#' package version 3.2.0. https://CRAN.R-project.org/package=maps
#'
#' Grosjean, Ph. (2016). svMisc: SciViews-R. UMONS, Mons, Belgium.
#' http://www.sciviews.org/SciViews-R.
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @author Ardern Hulme-Beaman
#' @export


IDbyDistanceDistInputCCVPar <- function(LatLongs, DistDataMat, Verbose=TRUE, ProvConfidence=0.95, Method=c('Spearman', 'Pearson')){

  UserInputAssessment(LatLongs, RefDistMat=DistDataMat, Method, RefData = NULL, DistVec = NULL)


  #making LatLongs a dataframe
  LatLongs <- as.data.frame(LatLongs)
  colnames(LatLongs) <- c("Lats", "Longs")


  GeoDist2PointPar <- function(RefLatLongs, TargetLatLong){
    GeographicDist <- apply(X = RefLatLongs[,2:1], MARGIN = 1, FUN = geosphere::distHaversine, p2 = TargetLatLong[2:1], r=6378.137)
    return(GeographicDist)
  }


  CorValCal <- function(ParRefLatLongs, ParTarLatLongs, ParDistDataMat, ParMethod){
    #ParRefLatLongs=LatLongs[-1,]; ParTarLatLongs=LatLongs[1,]; ParDistDataMat=DistDataMat[1,-1]; ParMethod=Method
    GeographicDistPar <- GeoDist2PointPar(RefLatLongs = ParRefLatLongs, TargetLatLong = ParTarLatLongs)

    if (ParMethod=='Spearman'){
      #running the correlation to generate r
      CorRes <- suppressWarnings(stats::cor.test(x = chr2nu(ParDistDataMat), y = GeographicDistPar, method = "spearman"))
    } else if (ParMethod=='Pearson'){
      CorRes <- suppressWarnings(stats::cor.test(x = chr2nu(ParDistDataMat), y = GeographicDistPar, method = "pearson"))
    }

    results <- c(ParTarLatLongs$Lats, ParTarLatLongs$Longs, CorRes$estimate)
    return(results)
  }

  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)

  b <- a <- 1
  CoordsHeat <- foreach::foreach(a = 1:dim(LatLongs)[1], b = -1:-dim(LatLongs)[1], .combine = rbind, .packages='GeoOrigins') %dopar% {
    CorValCal(ParRefLatLongs=LatLongs[b,], ParTarLatLongs=LatLongs[a,], ParDistDataMat=DistDataMat[a,b], ParMethod = Method)
  }

  parallel::stopCluster(clust)


  #converting the results into a data frame
  CoordsHeat <- as.data.frame(CoordsHeat, row.names = 1:dim(CoordsHeat)[1])

  #naming the variables
  names(CoordsHeat) <- c("Lats", "Longs", "Cor")

  ProvCor <- stats::quantile(CoordsHeat$Cor, 1-ProvConfidence)

  if (Verbose==TRUE){

    ProvResults <- list(Provenancing.Correlation=as.numeric(ProvCor), CCV.Cor.Vals=CoordsHeat$Cor)
    names(ProvResults)[1] <- paste(names(ProvResults)[1], ".", ProvConfidence*100, "%.Confidence", sep="")

    return(ProvResults)
  } else {
    ProvResults <- list(Provenancing.Correlation=as.numeric(ProvCor))
    names(ProvResults)[1] <- paste(names(ProvResults)[1], ".", ProvConfidence*100, "%.Confidence", sep="")

    return(ProvResults)
  }


}




#' Spatial Provenancing Correct Cross-Validation calculation from raw data using parallel processing
#'
#' This function takes raw variables of reference specimens with known spatial origins
#' and uses euclidean distances to calculate the correllation value that would be required to
#' correctly provenancing them in a test. This is achieved by a leave-one-out procedure. Shape variables
#' can be specified and if so Procrustes distances can be calculated. If a cross-validation procedure
#' that removes more than one specimen from the reference dataset is desired then it is recommended that
#' the validate function of the IDbyDistanceRawData method be used with a loop.
#' @param ProvConfidence is a value between 0 and 1 indicating the confidence level that is desired for spatial provenancing.
#' @inheritParams IDbyDistanceRawData
#' @return If Verbose is set to FALSE then a list of a single object containing the correlation value at the required confidence interval is returned. If Verbose is set to TRUE then a list is returned with two objects: the first is the correlation value at the required confidence interval; the second a dataframe of coordinates and the spatial-trait correlation values at the true locations of each specimen.
#' @details When used for shape data and for Procrustes distances this function makes use of the \code{\link[shapes]{procGPA}} and \code{\link[shapes]{procdist}} functions from the \code{shapes} package. When Euclidean distances are employed the \code{\link[stats]{dist}} function of the base \code{stats} package is used.
#' This method also makes use of the \code{\link[stats]{cor.test}} function from the \code{stats} package. When the \code{PrintProg} is set to TRUE, the \code{\link[svMisc]{progress}} function of the \code{svMisc} package is used.
#' The map plotting of this function makes use of the functions of the \code{maps} package.
#'
#' @section Citations:
#'
#' Original S code by Richard A. Becker, Allan R. Wilks. R version by Ray Brownrigg.
#' Enhancements by Thomas P Minka and Alex Deckmyn. (2017). maps: Draw Geographical Maps. R
#' package version 3.2.0. https://CRAN.R-project.org/package=maps
#'
#' Ian L. Dryden (2016). shapes: Statistical Shape Analysis. R package version 1.1-13.
#' https://CRAN.R-project.org/package=shapes
#'
#' Grosjean, Ph. (2016). svMisc: SciViews-R. UMONS, Mons, Belgium.
#' http://www.sciviews.org/SciViews-R.
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @author Ardern Hulme-Beaman
#' @export


IDbyDistanceRawDataCCVPar <- function(LatLongs, RefData, ShapeData=TRUE, ShapeDim=2, DistMethod=c("Euc", "Proc"), Verbose=TRUE, ProvConfidence=0.95, Method=c('Pearson', 'Spearman')){

  UserInputAssessment(LatLongs, RefData, Method, RefDistMat = NULL, DistVec = NULL)


  #making LatLongs a dataframe
  LatLongs <- as.data.frame(LatLongs)
  colnames(LatLongs) <- c("Lats", "Longs")



  #organising data for ease of analysis
  #combining ref and target specimens with target first
  if (ShapeData==TRUE){
    gpaRes <- shapes::procGPA(Mat2Array(RefData, LMdim = ShapeDim))
    TotalShape <- Array2Mat(gpaRes$rotated)
  } else {
    TotalShape <- RefData
  }

  #calculating euclidean distances between specimens
  #and then extracting distances to target specimen only
  if (DistMethod=="Euc"){
    ShapeDistMat <- as.matrix(stats::dist(TotalShape))
  } else if (DistMethod=="Proc" && ShapeData==TRUE){
    ShapeDistMat <- ProcDistanceTablePar(Mat2Array(TotalShape, LMdim=ShapeDim))
  } else if (DistMethod=="Proc" && ShapeData==FALSE){
    stop("Error: Procrustes distance selected, but ShapeData argument is set to FALSE")
  }


  GeoDist2PointPar <- function(RefLatLongs, TargetLatLong){
    GeographicDist <- apply(X = RefLatLongs[,2:1], MARGIN = 1, FUN = geosphere::distHaversine, p2 = TargetLatLong[2:1], r=6378.137)
    return(GeographicDist)
  }

  CorValCal <- function(ParRefLatLongs, ParTarLatLongs, ParDistDataMat, ParMethod){
    GeographicDistPar <- GeoDist2PointPar(RefLatLongs = ParRefLatLongs, TargetLatLong = ParTarLatLongs)
    #running the correlation to generate r

    if (ParMethod=='Spearman'){
      #running the correlation to generate r
      CorRes <- suppressWarnings(stats::cor.test(x = ParDistDataMat, y = GeographicDistPar, method = "spearman"))
    } else if (ParMethod=='Pearson'){
      CorRes <- suppressWarnings(stats::cor.test(x = ParDistDataMat, y = GeographicDistPar, method = "pearson"))
    }

    results <- c(ParTarLatLongs$Lats, ParTarLatLongs$Longs, CorRes$estimate)
    return(results)
  }

  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)

  b <- a <- 1
  CoordsHeat <- foreach::foreach(a = 1:dim(LatLongs)[1], b = -1:-dim(LatLongs)[1], .combine = rbind, .packages='GeoOrigins') %dopar% {
    CorValCal(ParRefLatLongs=LatLongs[b,], ParTarLatLongs=LatLongs[a,], ParDistDataMat=ShapeDistMat[a,b], ParMethod = Method)
  }

  parallel::stopCluster(clust)

  #converting the results into a data frame
  CoordsHeat <- as.data.frame(CoordsHeat, row.names = 1:dim(CoordsHeat)[1])



  #naming the variables
  names(CoordsHeat) <- c("Lats", "Longs", "Cor")

  ProvCor <- stats::quantile(CoordsHeat$Cor, 1-ProvConfidence)

  if (Verbose==TRUE){

    ProvResults <- list(Provenancing.Correlation=as.numeric(ProvCor), CCV.Cor.Vals=CoordsHeat$Cor)
    names(ProvResults)[1] <- paste(names(ProvResults)[1], ".", ProvConfidence*100, "%.Confidence", sep="")

    return(ProvResults)
  } else {
    ProvResults <- list(Provenancing.Correlation=as.numeric(ProvCor))
    names(ProvResults)[1] <- paste(names(ProvResults)[1], ".", ProvConfidence*100, "%.Confidence", sep="")

    return(ProvResults)
  }


}




#' Returns geographic boundaries identified from trait distances using parallel processing
#'
#' This function takes the distances among all reference specimens and uses a process similar to
#' leave-one-out correct cross-validation to identify likely spatial trait boundaries. This method
#' only allows for the inclusion of specimens with specific known locations (i.e. specific
#' latitude-longitude coordinates). This functions requires a distance input, which allows the user
#' to input any desired dissimilarity or distance metrics as appropriate for the data. This function
#' is a version of the \code{Boundaryfinder} function that makes use of the \code{foreach} package
#' @param RefDistMat is a square matrix of pairwise distances among all reference specimens
#' @param ExpandMap is a vector of 2 elements for expanding the plotting region of the map. The first element expands the latitudinal area and the second element expands the longitudinal area.
#' @param StartPoint is an iteger that denotes the specimen to start the process on. As the method needs to cycle through all known specimens in the database, this can take some time. If the process is stopped for whatever reason the process can be picked up again by adjusting the StartPoint to the specimen number that it had previously finished on.
#' @param RefIDs is a vector of the unique identifiers for each of the reference specimens in the reference dataset. These values will be used for naming the files in the datadump fololder and also for matching up with the returned summary results. Default is set to NULL and if this is not populated then reference data is worked through consecutively, naming the datadump files in consecutive order.
#' @param IgnorePrompts default is set to FALSE, but if set to TRUE queries such as those confirming the location of the datadump will be surpressed.
#' @inheritParams IDbyDistanceDistInput
#' @return An dataframe of all values for every sampled grid reference for every specimen is returned.
#' @details The resulting array can be transferred to the \code{PlotBoundaries} function for plotting.
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#'
#' @author Ardern Hulme-Beaman
#' @export


BoundaryFinderPar <- function(LatLongs, RefDistMat=matrix(), LongRange, LatRange, RangeSamp=10, ExpandMap=c(0,0), StartPoint=1, RefIDs=NULL, IgnorePrompts=FALSE, Method = c('Spearman', 'Pearson'), PacificCent=FALSE){

  #LatLongs = WolfM1Data$Lat.Long[1:10,]
  #RefDistMat = WolfM1Data$HeatDistMat[1:10,1:10]
  #LongRange = Long.Range
  #LatRange = Lat.Range
  #RangeSamp = LowResRsamp.Rp
  #ExpandMap = c(0,0)
  #RefIDs = WolfM1Data$info$ID
  #IgnorePrompts = TRUE
  #Method = 'Spearman'

  UserInputAssessment(LatLongs, RefDistMat, Method, RefData = NULL, DistVec = NULL)


  #making LatLongs a dataframe
  LatLongs <- as.data.frame(LatLongs)
  colnames(LatLongs) <- c("Lats", "Longs")


  #creating a range that will cover the whole geographic area of interest
  #this is for the function to loop through later on
  if (length(RangeSamp)==1){
    LongSamp <- RangeSamp
    LatSamp <- RangeSamp
  } else if (length(RangeSamp)==2){
    LongSamp <- RangeSamp[2]
    LatSamp <- RangeSamp[1]
  } else if (length(RangeSamp)>2){
    stop("too many dimensions in RangeSamp")
  }

  LongRangeSteps <- (LongRange[2]-LongRange[1])/(LongSamp-2)
  LatRangeSteps <- (LatRange[2]-LatRange[1])/(LatSamp-2)


  #this output for Lat/Long ways provides what the loop should sequence through
  Longways <- c(LongRange[1]-LongRangeSteps, seq(LongRange[1], LongRange[2], by = LongRangeSteps), LongRange[2]+LongRangeSteps)
  if(sum(Longways<c(-179))>0){
    Longways <- Longways[-which(Longways<c(-179))]
  }
  if(sum(Longways>c(179))>0){
    Longways <- Longways[-which(Longways>c(179))]
  }


  Latways <- c(LatRange[1]-LatRangeSteps, seq(LatRange[1], LatRange[2], by = LatRangeSteps), LatRange[2]+LatRangeSteps)
  if(sum(Latways<c(-90))>0){
    Latways <- Latways[-which(Latways<c(-90))]
  }
  if(sum(Latways>c(90))>0){
    Latways <- Latways[-which(Latways>c(90))]
  }


  if (is.null(RefIDs)){
    RefIDs <- 1:dim(RefDistMat)[1]
  }


  if (sum(rownames(RefDistMat)==c(1:length(rownames(RefDistMat))))==length(rownames(RefDistMat))){
    ArrayDimNames <- list(Latways, Longways, rownames(RefDistMat))
    if (IgnorePrompts==FALSE){
      warning("Distance matrix rownames appear match numeric order, so order number has been used as IDs for array dimensions.")
    }
  } else if (length(unique(rownames(RefDistMat)))==length(rownames(RefDistMat))){
    ArrayDimNames <- list(Latways, Longways, 1:length(rownames(RefDistMat)))
    if (IgnorePrompts==FALSE){
      print(x = "Distance matrix rownames appear to be unique and have been used as IDs for array dimensions.")
    }
  } else {
    ArrayDimNames <- list(Latways, Longways, 1:length(rownames(RefDistMat)))
    if (IgnorePrompts==FALSE){
      print(x = "Distance matrix rownames appear to have duplicates, so order number has been used as IDs for array dimensions.")
    }
  }

  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)



  IDbyDistanceDistInputPar <- function(LatLongsPar, DistDataVecPar, LongRangePar, LatRangePar, RangeSamp=10, ParMethod=Method, PacificCent=PacificCent){
    #LatLongsPar = LatLongs[-1,]; DistDataVecPar = RefDistMat[-1,1]; LongRangePar = LongRange; LatRangePar = LatRange; RangeSamp; ParMethod = Method

    #making LatLongs a dataframe
    LatLongsPar <- as.data.frame(LatLongsPar)
    colnames(LatLongsPar) <- c("Lats", "Longs")

    ShapeDist <- DistDataVecPar

    #creating an empty object to be populated by results
    CoordsHeat <- NULL


    #this function can come up with correlation values across the entire map
    #here we have the first option which is useful for the validation process
    #by doing this first option for all the specimens (in a loop) we can build a distribution
    #of the correlation values that will correctly cover the specimens true location
    #but see IDbyDistanceRawDataCCV function for looping to be done automatically

    #creating a range that will cover the whole geographic area of interest
    #this is for the function to loop through later on

    if (length(RangeSamp)==1){
      LongSamp <- RangeSamp
      LatSamp <- RangeSamp
    } else if (length(RangeSamp)==2){
      LongSamp <- RangeSamp[2]
      LatSamp <- RangeSamp[1]
    } else if (length(RangeSamp)>2){
      stop("too many dimensions in RangeSamp")
    }

    LongRangeSteps <- (LongRangePar[2]-LongRangePar[1])/(LongSamp-2)
    LatRangeSteps <- (LatRangePar[2]-LatRangePar[1])/(LatSamp-2)

    #this output for Lat/Long ways provides what the loop should sequence through
    if (PacificCent==TRUE){
      MidRange <- seq(LongRange[1], LongRange[2]+360, by = LongRangeSteps*-1)
      #MidRange[which(MidRange>=180)] <- MidRange[which(MidRange>=180)]-360
      Longways <- c(LongRange[1]+LongRangeSteps, MidRange, LongRange[2]+(LongRangeSteps*-1)+360)
    } else {
      Longways <- c(LongRange[1]-LongRangeSteps, seq(LongRange[1], LongRange[2], by = LongRangeSteps), LongRange[2]+LongRangeSteps)
    }


    Latways <- c(LatRangePar[1]-LatRangeSteps, seq(LatRangePar[1], LatRangePar[2], by = LatRangeSteps), LatRangePar[2]+LatRangeSteps)

    if(sum(Latways<c(-90))>0){
      Latways <- Latways[-which(Latways<c(-90))]
    }
    if(sum(Latways>c(90))>0){
      Latways <- Latways[-which(Latways>c(90))]
    }


    CorMatrixRes <- matrix(NA, nrow = length(Latways), ncol = length(Longways))
    rownames(CorMatrixRes) <- Latways
    colnames(CorMatrixRes) <- Longways

    GeoDist2PointPar <- function(RefLatLongs, TargetLatLong){
      GeographicDist <- apply(X = RefLatLongs[,2:1], MARGIN = 1, FUN = geosphere::distHaversine, p2 = TargetLatLong[2:1], r=6378.137)
      return(GeographicDist)
    }
    #carrying out iterative analyses across the geographic region defined by LongRange and LatRange
    for (i in 1:length(Longways)){
      for (j in 1:length(Latways)){
        #i <- 1
        #j <- 1

        coord <- c(Latways[j], Longways[i])

        GeographicDist <- GeoDist2PointPar(RefLatLongs = LatLongsPar, TargetLatLong = coord)

        if (ParMethod=='Spearman'){
          #running the correlation to generate r
          CorRes <- suppressWarnings(stats::cor.test(x = ShapeDist, y = GeographicDist, method = "spearman"))
        } else if (ParMethod=='Pearson'){
          CorRes <- suppressWarnings(stats::cor.test(x = ShapeDist, y = GeographicDist, method = "pearson"))
        }
        #populating results matrix
        CorMatrixRes[j,i] <- CorRes$estimate

      }


    }

    CoordsHeat <- CorMatrixRes
    return(CoordsHeat)
  }



  b <- a <- 1
  IdentificationRange <- foreach::foreach(a = StartPoint:dim(RefDistMat)[1], b = -StartPoint:-dim(RefDistMat)[1]) %dopar% {
    IDbyDistanceDistInputPar(LatLongsPar = LatLongs[b,], DistDataVecPar = RefDistMat[b,a], LongRangePar = LongRange, LatRangePar = LatRange, RangeSamp, ParMethod = Method, PacificCent = PacificCent)
  }

  parallel::stopCluster(clust)

  ProvArrayRes <- array(as.numeric(unlist(IdentificationRange)), dim=c(dim(IdentificationRange[[1]])[1], dim(IdentificationRange[[1]])[2], length(IdentificationRange)), dimnames = ArrayDimNames)

  return(list(RawCorData=ProvArrayRes))

}


#' Parallel processing calculation of trait boundary overlaps and statistics from Boundaryfinding R object
#'
#' This function takes the results of the \code{BoundaryFinder} function that have been
#' stored in an array and returns the statistics of the area coverage of boundaries
#' found for each specimen based on the provided r threshold. This function is particularly
#' useful when using the parallel processing function \code{BoundaryFinderPar}, which does not
#' return the Provenancing area values or percentages that are returned by the \code{BoundaryFInder}
#' function when the \code{verbose} aregument is set to \code{TRUE}. It is the parallel processing equivalent of
#' the \code{TraitBoundaryStats} function.
#' @inheritParams TraitBoundaryStats
#'
#'
#' @return This function returns a list of two objects: 1. the area of the convex hull of the distribution of samples and 2. a dataframe of two columns with the area returned of reach specimen in the identification process and the corresponding percentage that area represents of the sample distribution.
#'
#' @author Ardern Hulme-Beaman
#'
#' @export



TraitBoundaryStatsPar <- function(RawCorArray, LatLongs, RefIDs, PlotValCor){

  SpecimenLoc <- cbind(chr2nu(LatLongs$Longs), chr2nu(LatLongs$Lats))
  DistPol <-  sp::Polygon(SpecimenLoc[grDevices::chull(SpecimenLoc),])
  DistPols <-  sp::Polygons(list(DistPol),1)
  DistSpatPol <-  sp::SpatialPolygons(list(DistPols))
  sp::proj4string(DistSpatPol) <- "+proj=longlat +ellps=WGS84"
  DistSpatPol$area <- raster::area(DistSpatPol)


  Construct_contourPar <- function(LatLongs) {
    #LatLongs = as.data.frame(CoordsHeat)
    x <-  cbind(as.numeric(as.character(LatLongs$Longs)), as.numeric(as.character(LatLongs$Lats)))
    KSres <- ks::kde(x=x, H=ks::Hpi(x=x), compute.cont=TRUE)
    contour.95 <- with(KSres, grDevices::contourLines(x=eval.points[[1]],y=eval.points[[2]],z=estimate,levels=cont["5%"])[[1]])
    return(contour.95)
  }

  ParAreaFinder <- function(RawCorArrayPar, DistSpatPolPar, PlotValCorPar){
    #RawCorArrayPar = RawCorArray[,,1]; DistSpatPolPar = DistSpatPol; PlotValCorPar = PlotValCor

    CoordsHeat <- cbind(expand.grid(dimnames(RawCorArrayPar)[[1]], dimnames(RawCorArrayPar)[[2]]), c(RawCorArrayPar))

    colnames(CoordsHeat) <- c('Lats', 'Longs', 'Cor')

    Select.95.conf <- which(CoordsHeat$Cor>PlotValCor)

    if (length(Select.95.conf)!=0){

      ApproxOrigin <- CoordsHeat[Select.95.conf,]

      Latvar <- stats::var(as.numeric(as.character(ApproxOrigin$Lats)))
      Longvar <- stats::var(as.numeric(as.character(ApproxOrigin$Longs)))


      if (is.na(Latvar) || Latvar==0 || Longvar==0){
        EstDist <- NULL
      } else {
        contour.95 <-  Construct_contourPar(ApproxOrigin[,1:2])
        #polygon(contour.95$x, contour.95$y, col=transpar(Colour ='gold', alpha = 10), border='black', lwd=1)

        EstimatedPol <-  sp::Polygon(cbind(contour.95$x, contour.95$y))
        EstimatedPols <-  sp::Polygons(list(EstimatedPol),1)
        EstSpatPol <-  sp::SpatialPolygons(list(EstimatedPols))
        #plot(DistSpatPol)
        #plot(EstSpatPol, add=TRUE)
        sp::proj4string(EstSpatPol) <- "+proj=longlat +ellps=WGS84"


        EstDist <- raster::intersect(EstSpatPol, DistSpatPolPar)

      }

      if (is.null(EstDist)){
        NarrowedRanges <- c(0, 0)
      } else if (attr(class(EstDist), "package")!='sp'){
        NarrowedRanges <- c(0, 0)
      } else {
        EstDist$area <- raster::area(EstDist)
        NarrowedRanges <- c(EstDist$area/DistSpatPolPar$area, EstDist$area)
      }

    } else {
      NarrowedRanges <- c(1, DistSpatPolPar$area)
    }

    return(NarrowedRanges)
  }


  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)

  a <- 1
  AreaRes <- foreach::foreach(a = 1:dim(RawCorArray)[3], .combine = cbind, .packages=c('GeoOrigins', 'ks')) %dopar% {
    ParAreaFinder(RawCorArrayPar = RawCorArray[,,a], DistSpatPolPar = DistSpatPol, PlotValCorPar = PlotValCor)
  }

  parallel::stopCluster(clust)

  ResTableDimNames <- list(c('Percent.Overlap', 'Area.Overlap.in.m'), RefIDs)
  ResTable <- AreaRes
  dimnames(ResTable) <- ResTableDimNames

  return(ResTable)
}

