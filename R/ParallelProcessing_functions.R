


#' Spatially provenance a specimen from distance data adapted for use in parallel processing
#'
#' This function takes the distances from an unknown specimen to all reference specimens
#' and uses these distances to calculate a likely spatial provenance. Note that this
#' procedured can only be applied to one unknown specimen at a time. This distance input function
#' allows the user to input any dissimilarity or distance desired and as appropriate for the data.
#' The function has two applications either: calculating a specimens' provenance, or alternatively
#' it can be used to calculate the minimum correlation coefficient needed to correctly identify a
#' known specimen at its true collection location. The second application of this function can work
#' as a correct cross-validation process if looped, but see IDbyDistance.DistInput.CCV function which
#' does this automatically in a leave-one-out process. However, if a corss-validation process that
#' removed more than one specimen from the reference dataset at a time is required then it is advised
#' that this be applied using this function.
#' @param Dist.data.vec.Par is a vector of distances from each reference specimen to the specimen of interest (either a specimen with unknown provenance or another reference specimen that the user is interested in validating the provenance of)
#' @param Lat.Longs.Par is a
#' @inheritParams IDbyDistance.RawData
#' @return If verbose is TRUE then a dataframe of all values for every sampled grid reference is returned. If verbose is FALSE then only those grid references with the highest correlation values are returned.
#' @details This method also makes use of the \code{cor.test} function from the \code{stats} package. When the \code{print.prog} is set to TRUE, the \code{progress} function of the \code{svMisc} package is used.
#' The map plotting of this function makes use of the functions of the \code{maps} package.
#'
#' @keywords Spatial provenancing
#' @keywords Spatial identification
#' @keywords internal
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#' @author Ardern Hulme-Beaman




IDbyDistance.DistInput.Par <- function(Lat.Longs.Par, Dist.data.vec.Par, LongRange.Par, LatRange.Par, Range.Samp=10, Validate=FALSE){
  #Lat.Longs = cbind(RatLatLongData$Lat, RatLatLongData$Long)[-1,]; Dist.data.vec = RatDistMat[-1,1]; LongRange = Long.Range; LatRange = Lat.Range; Range.Samp = R.Samp; verbose = FALSE; Validate = FALSE
  #Lat.Longs = NEFtInfo[-1,2:3]; Dist.data.vec = NEDisMat[-1,1]; LongRange = NEFt.Long.Range; LatRange = NEFt.Lat.Range; Range.Samp = NEFtR.Samp; verbose = FALSE; Validate = FALSE; plot.Val.cor = NEFtrThres$`Provenancing.Correlation.95%.Confidence`; plot.Prov = TRUE
  #Lat.Longs = RpraetorLatLong[-1,]; Dist.data.vec = RatDistMat[-1,1]; LongRange = Long.Range; LatRange = Lat.Range; Range.Samp = R.Samp; verbose = FALSE; Validate = FALSE; plot.Val.cor = rThres$`Provenancing.Correlation.95%.Confidence`; plot.Prov = TRUE

  #making Lat.Longs a dataframe
  Lat.Longs <- as.data.frame(Lat.Longs.Par)
  colnames(Lat.Longs) <- c("Lats", "Longs")


  Shape.Dist <- Dist.data.vec.Par


  #creating an empty object to be populated by results
  CoordsHeat <- NULL


  #this function can come up with correlation values across the entire map
  #here we have the first option which is useful for the validation process
  #by doing this first option for all the specimens (in a loop) we can build a distribution
  #of the correlation values that will correctly cover the specimens true location
  #but see IDbyDistance.RawData.ccv function for looping to be done automatically

  #creating a range that will cover the whole geographic area of interest
  #this is for the function to loop through later on

  if (length(Range.Samp)==1){
    Long.Samp <- Range.Samp
    Lat.Samp <- Range.Samp
  } else if (length(Range.Samp)==2){
    Long.Samp <- Range.Samp[2]
    Lat.Samp <- Range.Samp[1]
  } else if (length(Range.Samp)>2){
    stop("too many dimensions in Range.Samp")
  }

  Long.Range.Steps <- (LongRange.Par[2]-LongRange.Par[1])/(Long.Samp-2)
  Lat.Range.Steps <- (LatRange.Par[2]-LatRange.Par[1])/(Lat.Samp-2)

  #this output for Lat/Long ways provides what the loop should sequence through
  Longways <- c(LongRange.Par[1]-Long.Range.Steps, seq(LongRange.Par[1], LongRange.Par[2], by = Long.Range.Steps), LongRange.Par[2]+Long.Range.Steps)
  Latways <- c(LatRange.Par[1]-Lat.Range.Steps, seq(LatRange.Par[1], LatRange.Par[2], by = Lat.Range.Steps), LatRange.Par[2]+Lat.Range.Steps)

  cor.matrix.Res <- matrix(NA, nrow = length(Latways), ncol = length(Longways))
  rownames(cor.matrix.Res) <- Latways
  colnames(cor.matrix.Res) <- Longways

  Geo.Dist2Point.Par <- function(RefLatLongs, TargetLatLong){
    Geographic.Dists <- apply(X = RefLatLongs[,2:1], MARGIN = 1, FUN = geosphere::distHaversine, p2 = TargetLatLong[2:1])
    return(Geographic.Dists)
  }
  #carrying out iterative analyses across the geographic region defined by LongRange and LatRange
  for (i in 1:length(Longways)){
    for (j in 1:length(Latways)){
      #i <- Longways[1]
      #j <- Latways[1]

      coord <- c(Latways[j], Longways[i])

      Geographic.Dists <- Geo.Dist2Point.Par(RefLatLongs = Lat.Longs, TargetLatLong = coord)

      #running the correlation to generate r
      CorRes <- stats::cor.test(x = Shape.Dist, y = Geographic.Dists)

      #populating results matrix
      cor.matrix.Res[j,i] <- CorRes$estimate

    }


  }

  CoordsHeat <- cor.matrix.Res

  #if there is not a validating the data then we this means
  #we either have the validation result from a previous analyses and we can plot it
  #or we don't and we're not yet interested in it
  return(CoordsHeat)
}





#' Spatial Provenancing Correct Cross-Validation calculation from distance data using parallel processing
#'
#' This function takes pairwise distances among all reference specimens with known spatial origins
#' and uses those distances to calculate the correllation value that would be required to
#' correctly provenancing them in a test. This is achieved by a leave-one-out procedure. If a
#' cross-validation procedure that removes more than one specimen from the reference dataset is
#' desired then it is recommended that the validate function of the IDbyDistance.DistInput method be used
#' with the validate argument in a loop.
#' @param Dist.data.mat is a square matrix of pairwise distances among all reference specimens
#' @inheritParams IDbyDistance.RawData.CCV
#' @return If verbose is set to FALSE then a list of a single object containing the correlation value at the required confidence interval is returned. If verbose is set to TRUE then a list is returned with two objects: the first is the correlation value at the required confidence interval; the second a dataframe of coordinates and the spatial-trait correlation values at the true locations of each specimen.
#' @details This method also makes use of the \code{cor.test} function from the \code{stats} package. When the \code{print.prog} is set to TRUE, the \code{progress} function of the \code{svMisc} package is used.
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


IDbyDistance.DistInput.CCV.Par <- function(Lat.Longs, Dist.data.mat, verbose=TRUE, Prov.Confidence=0.95){
  #Dist.data.mat = RatDistMat

  #making Lat.Longs a dataframe
  Lat.Longs <- as.data.frame(Lat.Longs)
  colnames(Lat.Longs) <- c("Lats", "Longs")


  Geo.Dist2Point.Par <- function(RefLatLongs, TargetLatLong){
    Geographic.Dists <- apply(X = RefLatLongs[,2:1], MARGIN = 1, FUN = geosphere::distHaversine, p2 = TargetLatLong[2:1])
    return(Geographic.Dists)
  }


  CorValCal <- function(Par.Ref.Lat.Longs, Par.Tar.Lat.Longs, Par.Dist.data.mat){
    Geographic.Dists <- Geo.Dist2Point.Par(RefLatLongs = Par.Ref.Lat.Longs, TargetLatLong = Par.Tar.Lat.Longs)
    #running the correlation to generate r
    CorRes <- stats::cor.test(x = Par.Dist.data.mat, y = Geographic.Dists)
    results <- c(Par.Tar.Lat.Longs$Lats, Par.Tar.Lat.Longs$Longs, CorRes$estimate)
    return(results)
  }

  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)

  b <- a <- 1
  CoordsHeat <- foreach::foreach(a = 1:dim(Lat.Longs)[1], b = -1:-dim(Lat.Longs)[1], .combine = rbind) %dopar% {
    CorValCal(Par.Ref.Lat.Longs=Lat.Longs[b,], Par.Tar.Lat.Longs=Lat.Longs[a,], Par.Dist.data.mat=Dist.data.mat[a,b])
  }

  parallel::stopCluster(clust)


  #converting the results into a data frame
  CoordsHeat <- as.data.frame(CoordsHeat, row.names = 1:dim(CoordsHeat)[1])

  #naming the variables
  names(CoordsHeat) <- c("Lats", "Longs", "Cor")

  ProvCor <- stats::quantile(CoordsHeat$Cor, 1-Prov.Confidence)

  if (verbose==TRUE){

    ProvResults <- list(Provenancing.Correlation=as.numeric(ProvCor), CCV.Cor.Vals=CoordsHeat$Cor)
    names(ProvResults)[1] <- paste(names(ProvResults)[1], ".", Prov.Confidence*100, "%.Confidence", sep="")

    return(ProvResults)
  } else {
    ProvResults <- list(Provenancing.Correlation=as.numeric(ProvCor))
    names(ProvResults)[1] <- paste(names(ProvResults)[1], ".", Prov.Confidence*100, "%.Confidence", sep="")

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
#' the validate function of the IDbyDistance.RawData method be used with a loop.
#' @param Prov.Confidence is a value between 0 and 1 indicating the confidence level that is desired for spatial provenancing.
#' @inheritParams IDbyDistance.RawData
#' @return If verbose is set to FALSE then a list of a single object containing the correlation value at the required confidence interval is returned. If verbose is set to TRUE then a list is returned with two objects: the first is the correlation value at the required confidence interval; the second a dataframe of coordinates and the spatial-trait correlation values at the true locations of each specimen.
#' @details When used for shape data and for Procrustes distances this function makes use of the \code{procGPA} and \code{procdist} functions from the \code{shapes} package. When Euclidean distances are employed the \code{dist} function of the base \code{stats} package is used.
#' This method also makes use of the \code{cor.test} function from the \code{stats} package. When the \code{print.prog} is set to TRUE, the \code{progress} function of the \code{svMisc} package is used.
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


IDbyDistance.RawData.CCV.Par <- function(Lat.Longs, Ref.data, Shape.Data=TRUE, Shape.dim=2, DistMethod=c("Euc", "Proc"), verbose=TRUE, Prov.Confidence=0.95){

  #making Lat.Longs a dataframe
  Lat.Longs <- as.data.frame(Lat.Longs)
  colnames(Lat.Longs) <- c("Lats", "Longs")



  #organising data for ease of analysis
  #combining ref and target specimens with target first
  if (Shape.Data==TRUE){
    gpaRes <- shapes::procGPA(Mat2Array(Ref.data, LMdim = Shape.dim))
    total.shape <- Array2Mat(gpaRes$rotated)
  } else {
    total.shape <- Ref.data
  }

  #calculating euclidean distances between specimens
  #and then extracting distances to target specimen only
  if (DistMethod=="Euc"){
    Shape.Dist.Mat <- as.matrix(stats::dist(total.shape))
  } else if (DistMethod=="Proc" && Shape.Data==TRUE){
    Shape.Dist.Mat <- ProcDistanceTable.Par(Mat2Array(total.shape, LMdim=Shape.dim))
  } else if (DistMethod=="Proc" && Shape.Data==FALSE){
    stop("Error: Procrustes distance selected, but Shape.Data argument is set to FALSE")
  }


  Geo.Dist2Point.Par <- function(RefLatLongs, TargetLatLong){
    Geographic.Dists <- apply(X = RefLatLongs[,2:1], MARGIN = 1, FUN = geosphere::distHaversine, p2 = TargetLatLong[2:1])
    return(Geographic.Dists)
  }

  CorValCal <- function(Par.Ref.Lat.Longs, Par.Tar.Lat.Longs, Par.Dist.data.mat){
    Geographic.Dists <- Geo.Dist2Point.Par(RefLatLongs = Par.Ref.Lat.Longs, TargetLatLong = Par.Tar.Lat.Longs)
    #running the correlation to generate r
    CorRes <- stats::cor.test(x = Par.Dist.data.mat, y = Geographic.Dists)
    results <- c(Par.Tar.Lat.Longs$Lats, Par.Tar.Lat.Longs$Longs, CorRes$estimate)
    return(results)
  }

  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)

  b <- a <- 1
  CoordsHeat <- foreach::foreach(a = 1:dim(Lat.Longs)[1], b = -1:-dim(Lat.Longs)[1], .combine = rbind) %dopar% {
    CorValCal(Par.Ref.Lat.Longs=Lat.Longs[b,], Par.Tar.Lat.Longs=Lat.Longs[a,], Par.Dist.data.mat=Shape.Dist.Mat[a,b])
  }

  parallel::stopCluster(clust)

  #converting the results into a data frame
  CoordsHeat <- as.data.frame(CoordsHeat, row.names = 1:dim(CoordsHeat)[1])



  #naming the variables
  names(CoordsHeat) <- c("Lats", "Longs", "Cor")

  ProvCor <- stats::quantile(CoordsHeat$Cor, 1-Prov.Confidence)

  if (verbose==TRUE){

    ProvResults <- list(Provenancing.Correlation=as.numeric(ProvCor), CCV.Cor.Vals=CoordsHeat$Cor)
    names(ProvResults)[1] <- paste(names(ProvResults)[1], ".", Prov.Confidence*100, "%.Confidence", sep="")

    return(ProvResults)
  } else {
    ProvResults <- list(Provenancing.Correlation=as.numeric(ProvCor))
    names(ProvResults)[1] <- paste(names(ProvResults)[1], ".", Prov.Confidence*100, "%.Confidence", sep="")

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
#' @param Ref.Dist.mat is a square matrix of pairwise distances among all reference specimens
#' @param Expand.map is a vector of 2 elements for expanding the plotting region of the map. The first element expands the latitudinal area and the second element expands the longitudinal area.
#' @param startpoint is an iteger that denotes the specimen to start the process on. As the method needs to cycle through all known specimens in the database, this can take some time. If the process is stopped for whatever reason the process can be picked up again by adjusting the startpoint to the specimen number that it had previously finished on.
#' @param Ref.IDs is a vector of the unique identifiers for each of the reference specimens in the reference dataset. These values will be used for naming the files in the datadump fololder and also for matching up with the returned summary results. Default is set to NULL and if this is not populated then reference data is worked through consecutively, naming the datadump files in consecutive order.
#' @param ignore.prompts default is set to FALSE, but if set to TRUE queries such as those confirming the location of the datadump will be surpressed.
#' @inheritParams IDbyDistance.DistInput
#' @return An dataframe of all values for every sampled grid reference for every specimen is returned.
#' @details The resulting array can be transferred to the \code{PlotBoundaries} function for plotting.
#'
#' @import doParallel
#' @import parallel
#' @import foreach
#'
#' @keywords Spatial boundary
#' @keywords Spatial cluster
#' @author Ardern Hulme-Beaman
#' @export


BoundaryFinder.Par <- function(Lat.Longs, Ref.Dist.mat=matrix(), LongRange, LatRange, Range.Samp=10, Expand.map=c(0,0), startpoint=1, Ref.IDs=NULL, ignore.prompts=FALSE){
  #Lat.Longs = Rpraetor$Lat.Long; Ref.Dist.mat = RatDistMat;LongRange = range(Rpraetor$Lat.Long$Long); LatRange = range(Rpraetor$Lat.Long$Lat); Range.Samp = c(20,10); plot.Val.cor = rThres$`Provenancing.Correlation.95%.Confidence`; Expand.map = c(0,0); Ref.IDs = 1:dim(Rpraetor$LMs)[3]; DataDump = FALSE
  #Expand.map=c(0,0); DataDump=TRUE; DataDumpPath=NA; startpoint=1; Ref.IDs=NULL; ignore.prompts=FALSE



  #making Lat.Longs a dataframe
  Lat.Longs <- as.data.frame(Lat.Longs)
  colnames(Lat.Longs) <- c("Lats", "Longs")


  #creating a range that will cover the whole geographic area of interest
  #this is for the function to loop through later on
  if (length(Range.Samp)==1){
    Long.Samp <- Range.Samp
    Lat.Samp <- Range.Samp
  } else if (length(Range.Samp)==2){
    Long.Samp <- Range.Samp[2]
    Lat.Samp <- Range.Samp[1]
  } else if (length(Range.Samp)>2){
    stop("too many dimensions in Range.Samp")
  }

  Long.Range.Steps <- (LongRange[2]-LongRange[1])/(Long.Samp-2)
  Lat.Range.Steps <- (LatRange[2]-LatRange[1])/(Lat.Samp-2)


  #this output for Lat/Long ways provides what the loop should sequence through
  Longways <- c(LongRange[1]-Long.Range.Steps, seq(LongRange[1], LongRange[2], by = Long.Range.Steps), LongRange[2]+Long.Range.Steps)
  Latways <- c(LatRange[1]-Lat.Range.Steps, seq(LatRange[1], LatRange[2], by = Lat.Range.Steps), LatRange[2]+Lat.Range.Steps)


  if (is.null(Ref.IDs)){
    Ref.IDs <- 1:dim(Ref.Dist.mat)[1]
  }

  if (sum(rownames(Ref.Dist.mat)==c(1:length(rownames(Ref.Dist.mat))))==length(rownames(Ref.Dist.mat))){
    ArrayDimNames <- list(Latways, Longways, rownames(Ref.Dist.mat))
    if (ignore.prompts==FALSE){
      readline(prompt = "\n Distance matrix rownames appear match numeric order, \n so order number has been used as IDs for DataDump file names or array dimensions. \n \n Press enter to continue.")
    }
  } else if (length(unique(rownames(Ref.Dist.mat)))==length(rownames(Ref.Dist.mat))){
    ArrayDimNames <- list(Latways, Longways, 1:length(rownames(Ref.Dist.mat)))
    if (ignore.prompts==FALSE){
      readline(prompt = "\n Distance matrix rownames appear to be unique \n and have been used as IDs for DataDump file names or array dimensions. \n \n Press enter to continue.")
    }
  } else {
    ArrayDimNames <- list(Latways, Longways, 1:length(rownames(Ref.Dist.mat)))
    if (ignore.prompts==FALSE){
      readline(prompt = "\n Distance matrix rownames appear to have duplicates, \n so order number has been used as IDs for DataDump file names or array dimensions. \n \n Press enter to continue.")
    }
  }

  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)

  b <- a <- 1
  IdentificationRange <- foreach::foreach(a = startpoint:dim(Ref.Dist.mat)[1], b = -startpoint:-dim(Ref.Dist.mat)[1]) %dopar% {
    IDbyDistance.DistInput.Par(Lat.Longs.Par = Lat.Longs[b,], Dist.data.vec.Par = Ref.Dist.mat[b,a], LongRange.Par = LongRange, LatRange.Par = LatRange, Range.Samp)
  }

  parallel::stopCluster(clust)

  Prov.Array.Res <- array(as.numeric(unlist(IdentificationRange)), dim=c(dim(IdentificationRange[[1]])[1], dim(IdentificationRange[[1]])[2], length(IdentificationRange)), dimnames = ArrayDimNames)

  return(list(RawCorData=Prov.Array.Res))

}






