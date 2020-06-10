

#' Spatially provenance a specimen from raw trait data
#'
#' This function takes the raw variables of an unknown specimen and reference specimens
#' and uses euclidean distances to calculate a likely spatial provenance. Note that this
#' procedure can only be applied to one unknown specimen at a time. Shape variables
#' can be specified and if so Procrustes distances can be calculated. It has two applications
#' either: calculating a specimens' provenance, or alternatively it can be used to calculate the
#' minimum correlation coefficient needed to correctly identify a known specimen at its true
#' collection location. The second application of this function can work as a correct cross-
#' validation process if looped, but see \code{IDbyDistanceRawDataCCV} function which does this automatically
#' in a leave-one-out process. However, if a cross-validation process that removed more than one specimen
#' from the reference dataset at a time is required then it is advised that this be applied using
#' this function.
#' @param LatLongs a matrix of n rows by 2 columns where n is the number of reference specimens in your dataset and the columns are Latitude and Longitude values in that order. These latitude-longitude coordinates should be of the locations of the reference specimens.
#' @param TargetData is a vector of unknown specimen data. If it is geometric morphometric data is should be a vector of superimposed coordinated in the format X1, Y1, X2, Y2... etc. (or a vector of other standardised variables that can appropriately have euclidean distances calculated between it and reference variables). NB if applied to a reference collection specimen be sure to remove it from the RefData dataset.
#' @param RefData is a matrix of Reference specimen data where the rows are the individual reference specimens and the columns are the variables in the same order as the TargetData vector.
#' @param ShapeData logical indicating whether the data is geometric morphometric shape data that requires superimposition. Default set to TRUE
#' @param ShapeDim integer either 2 or 3 to indicate the dimensions of landmark coordinates if the data is geometric morphometric data.
#' @param DistMethod determines what kind of distance calculation should be used, either Euclidean "Euc" or Procrustes "Proc". If the user wishes to use another distance or dissimilarity please use the IDbyDistanceDistInput function.
#' @param LongRange is a vector of 2 elements defining the maximum and minimum Longitude values that the provenancing method should explore. This will also define the mapping range in the final plotted output.
#' @param LatRange is a vector of 2 elements defining the maximum and minimum Latitude values that the provenancing method should explore. This will also define the mapping range in the final plotted output.
#' @param RangeSamp is an integer vector of 1 or 2 elements that defines the resolution of spatial sampling. If one element is provided then both the latitude and longitude ranges are equally and evenly sampled using this value. If 2 elements are provided they should be in the order of latitude longitude and each range will be evenly sampled with its respective value.
#' @param Verbose logical whether or not a matrix of spatial correlation values is returned or not. Default is set to TRUE.
#' @param PrintProg logical whether or not to print a progress bar. Default set to FALSE.
#' @param Validate logical whether or not to run a correct cross-validation analysis to find the lowest required correlation value for correct identification.
#' @param ValidLatLongs if the process is carried out on a specimen of known location `(e.g. Validate=TRUE)`, then the latitude longitude coordinates for that location should be provided here in that order.
#' @param PlotRes logical whether or not to plot the provenancing map with heat values of most likely spatial origin.
#' @param HeatHue numeric vector of 2 elements each between 0 and 1. The first should be the hue value on the HSV scale; the second value should be the level of transparency of the colour used.
#' @param TileSize numeric to dictate the pixel size of the heat mapping colour values.
#' @param PlotProv logical if the map should be printed with a polygon demarking a contour at a user defined correlation value.
#' @param PlotValCor numeric correlation value that is used to determine the most likely origin of the specimen. This value can be calculated by using the correct cross-validation method and identifying the correlation value that will correctly identify a desired percentage of specimens (e.g. 95\%).
#' @param Method determines what kind of correlation coefficient should be used, either "Spearman" or "Pearson". Spearman's ranked correlation coefficient does not assume a linear relationship between geographic and trait distances, whereas Pearson's coefficient does.
#' @return If Verbose is TRUE then a dataframe of all values for every sampled grid reference is returned. If Verbose is FALSE then only those grid references with the highest correlation values are returned.
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
#' @author Ardern Hulme-Beaman
#'
#' @examples
#' Range.Exp <- .5
#'
#' Long.Range <- c(floor(min(Rpraetor$Lat.Long$Long))
#'               -Range.Exp,ceiling(max(Rpraetor$Lat.Long$Long)+Range.Exp))
#' Lat.Range <- c(floor(min(Rpraetor$Lat.Long$Lat))
#'               -Range.Exp,ceiling(max(Rpraetor$Lat.Long$Lat)+Range.Exp))
#'
#' RpraetorDataMat <- Array2Mat(Rpraetor$LMs)
#'
#' rThres <- IDbyDistanceRawDataCCV(LatLongs = Rpraetor$Lat.Long,
#'               RefData = RpraetorDataMat,
#'               ShapeData=TRUE,
#'               ShapeDim=2,
#'               DistMethod= "Proc",
#'               Verbose = TRUE,
#'               ProvConfidence = .95,
#'               PrintProg = FALSE,
#'               Method = 'Spearman')
#'
#' R.Samp <- c(12, 42)
#'
#' IDbyDistanceRawData(LatLongs = Rpraetor$Lat.Long[-1,],
#'                  TargetData = RpraetorDataMat[1,],
#'                  RefData = RpraetorDataMat[-1,],
#'                  ShapeData = TRUE,
#'                  ShapeDim = 2,
#'                  DistMethod = "Proc",
#'                  LongRange = Long.Range,
#'                  LatRange = Lat.Range,
#'                  RangeSamp = R.Samp,
#'                  Verbose = FALSE,
#'                  Validate = FALSE,
#'                  PlotValCor = rThres$`Provenancing.Correlation.95%.Confidence`,
#'                  PlotProv = TRUE,
#'                  Method = 'Spearman')
#'
#' points(x = Rpraetor$Lat.Long$Long[1], y=Rpraetor$Lat.Long$Lat[1], col='blue', pch=16)
#'
#' @export
#'
#'

IDbyDistanceRawData <- function(LatLongs, TargetData, RefData, ShapeData=TRUE, ShapeDim=2, DistMethod=c("Euc", "Proc"), LongRange=c(0,0), LatRange=c(0,0), RangeSamp=10, Verbose=TRUE, PrintProg=FALSE, Validate= FALSE, ValidLatLongs, PlotRes=TRUE, HeatHue= c(.15, 1), TileSize=2, PlotProv=FALSE, PlotValCor, Method=c('Spearman', 'Pearson')){

  UserInputAssessment(LatLongs, RefData, Method, RefDistMat = 'skip', DistVec = 'skip')


  #making LatLongs a dataframe
  LatLongs <- as.data.frame(LatLongs)
  colnames(LatLongs) <- c("Lats", "Longs")



  #organising data for ease of analysis
  #combining ref and target specimens with target first
  if (ShapeData==TRUE){
    TotalShape.raw <- rbind(TargetData, RefData)
    gpaRes <- shapes::procGPA(Mat2Array(TotalShape.raw, LMdim = ShapeDim))
    TotalShape <- Array2Mat(gpaRes$rotated)
  } else {
    TotalShape <- rbind(TargetData, RefData)
  }



  #calculating euclidean distances between specimens
  #and then extracting distances to target specimen only
  if (DistMethod=="Euc"){
    ShapeDist <- as.matrix(stats::dist(TotalShape))[-1,1]
  } else if (DistMethod=="Proc" && ShapeData==TRUE){
    ShapeDist <- apply(X = Mat2Array(TotalShape[-1,], LMdim=ShapeDim), MARGIN = 3, FUN = shapes::procdist, y=matrix(TotalShape[1,], nrow = length(TotalShape[1,])/ShapeDim, ncol = ShapeDim, byrow = TRUE))
  } else if (DistMethod=="Proc" && ShapeData==FALSE){
    stop("Error: Procrustes distance selected, but ShapeData argument is set to FALSE")
  }






  #creating an empty object to be populated by results
  CoordsHeat <- NULL

  #this function has two operations
  #either it can come up with the correlation value at a specific point (Validate=TRUE)
  #or it can come up with correlation values across the entire map
  #here we have the first option which is useful for the validation process
  #by doing this first option for all the specimens (in a loop) we can build a distribution
  #of the correlation values that will correctly cover the specimens true location
  #but see IDbyDistanceRawDataCCV function for looping to be done automatically
  if (Validate==TRUE){

    names(ValidLatLongs) <- c("Lats", "Longs")

    #This generates all the distances from the point on the map to all the specimen locations
    GeographicDist <- GeoDist2Point(RefLatLongs = LatLongs, TargetLatLong = c(ValidLatLongs$Lats, ValidLatLongs$Longs))

    #just checking the correlation visually
    #plot((ShapeDist), (GeographicDist))


    if (Method=='Spearman'){
      #running the correlation to generate r
      CorRes <- suppressWarnings(stats::cor.test(x = ShapeDist, y = GeographicDist, method = "spearman"))
    } else if (Method=='Pearson'){
      CorRes <- suppressWarnings(stats::cor.test(x = ShapeDist, y = GeographicDist, method = "pearson"))
    }

    #combining the individual results and organising them
    results <- c(ValidLatLongs, CorRes$estimate)

    #adding the results of this round to the previous results
    CoordsHeat <- rbind(CoordsHeat, results)
  } else if (Validate==FALSE){
    #this is the second of the above options and not for validation
    #here the function goes through the entire map and generated correlation values
    #the level of detail is set by "RangeSamp"
    #for example if RangeSamp is set to 1 the following procedure
    #will generate correlation values at every degree of the map
    #if the RangeSamp is set to 2 it will generate for every other degree
    #this can speed things up for quick and dirty tests

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
    Latways <- c(LatRange[1]-LatRangeSteps, seq(LatRange[1], LatRange[2], by = LatRangeSteps), LatRange[2]+LatRangeSteps)


    #carrying out iterative analyses across the geographic region defined by LongRange and LatRange
    for (i in Longways){
      for (j in Latways){
        coord <- c(j, i)

        GeographicDist <- GeoDist2Point(RefLatLongs = LatLongs, TargetLatLong = coord)

        if (Method=='Spearman'){
          #running the correlation to generate r
          CorRes <- suppressWarnings(stats::cor.test(x = ShapeDist, y = GeographicDist, method = "spearman"))
        } else if (Method=='Pearson'){
          CorRes <- suppressWarnings(stats::cor.test(x = ShapeDist, y = GeographicDist, method = "pearson"))
        }

        #combining the individual results and organising them
        results <- c(as.character(j),i, CorRes$estimate)

        #adding the results of this round to the previous results
        CoordsHeat <- rbind(CoordsHeat, results)
      }

      if (PrintProg==TRUE){
        svMisc::progress(value = which(Longways==i), max.value = length(Longways), progress.bar = FALSE)
        Sys.sleep(0.01)

      }

    }
  }


  #converting the results into a data frame
  CoordsHeat <- as.data.frame(CoordsHeat, row.names = 1:dim(CoordsHeat)[1])

  #naming the variables
  names(CoordsHeat) <- c("Lats", "Longs", "Cor")


  #if there is not a validating the data then we this means
  #we either have the validation result from a previous analyses and we can plot it
  #or we don't and we're not yet interested in it
  if (Validate==FALSE){
    if (PlotRes==TRUE){
      maps::map("world", xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="black", bg=graphics::par(bg="white"))


      #creating colour scale from max and min correlation based on which variable we're using
      CoordsHeatNum <- chr2nu(CoordsHeat$Cor)
      OriginLoc <- CoordsHeat[which(CoordsHeatNum==max(CoordsHeatNum)),]


      CoordsHeatscaled <- (CoordsHeatNum-min(CoordsHeatNum))/(max(CoordsHeatNum)-min(CoordsHeatNum))

      CoordsHeats <- grDevices::hsv(h = HeatHue[1], v = 1, s = CoordsHeatscaled, alpha = HeatHue[2])

      #plotting the correlations
      graphics::points(x = as.character(CoordsHeat$Long), y = as.character(CoordsHeat$Lat), pch=15, col=CoordsHeats,  cex=TileSize)

      #here if we want to plot a polygon of the region that
      #approximates the region the specimens came from (with whatever level of confidence we have selected)
      if (PlotProv==TRUE){
        Select.95.conf <- which(chr2nu(CoordsHeat$Cor)>PlotValCor)

        ApproxOrigin <- CoordsHeat[Select.95.conf,]

        Latvar <- stats::var(chr2nu(ApproxOrigin$Lats))
        Longvar <- stats::var(chr2nu(ApproxOrigin$Longs))

        if (is.na(Latvar)){
          RthresMessage <- paste("A provenancing region was not identified at this r threshold. R threshold set to:", PlotValCor, sep = "   ")
          MaxCorMessage <- paste("The maximum correlation value acheived was:", max(chr2nu(CoordsHeat$Cor)), sep = "   ")
          WarningMessage <- "This may be because the resolution used is too low so the likely origin region has been overlooked or alternatively the specimen could not be successfully identified because the reference material does not suffieciently reflect the morphology of the unknown specimen. Please adjust the sampling resolution by changing the R.Samp argument or change the r threshold or consider the specimen unidentifiable to a given region."
          warning(paste(RthresMessage, MaxCorMessage, WarningMessage, sep = " "))

        } else if (Latvar==0 || Longvar==0){
          warning("The resolution used for identifying a region of identification is too low to plot a polygon of the likely region of origin. Therefore, the grid squares that were identified by the r threshold as a possible region of origin have been highlighted. Please set the R.Samp argument to a higher value if a polygon of the most likely region of origin is desired.")
          polycol <- grDevices::hsv(h = HeatHue[1], s = 1, v = .8, alpha = HeatHue[2])
          graphics::points(x = as.character(ApproxOrigin$Longs), y = as.character(ApproxOrigin$Lats), pch=22, bg=polycol,  cex=TileSize+.1)
        } else {

          contour.95 <-  Construct_contour(ApproxOrigin)
          polycol <- grDevices::hsv(h = HeatHue[1], s = 1, v = .8, alpha = HeatHue[2])

          graphics::polygon(contour.95$x, contour.95$y, col=NA, border=polycol, lwd=2)

        }


      }

      maps::map("world", xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="black", bg=graphics::par(bg="white"), add=T)

      graphics::points(x = as.character(LatLongs$Long), y = as.character(LatLongs$Lat), pch=23, bg='orange',  cex=1)
      maps::map.axes()
    }
  }


  if (Verbose==TRUE|Validate==TRUE){
    return(CoordsHeat)
  } else {
    OriginLocCor <- CoordsHeat[which(CoordsHeat$Cor==max(chr2nu(CoordsHeat$Cor), na.rm = TRUE)),]

    return(list(Cor=OriginLocCor))
  }


}



#' Spatial Provenancing Correct Cross-Validation calculation from raw data
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
#' @author Ardern Hulme-Beaman
#'
#' @examples
#' Range.Exp <- .5
#'
#' Long.Range <- c(floor(min(Rpraetor$Lat.Long$Long))
#'               -Range.Exp,ceiling(max(Rpraetor$Lat.Long$Long)+Range.Exp))
#' Lat.Range <- c(floor(min(Rpraetor$Lat.Long$Lat))
#'               -Range.Exp,ceiling(max(Rpraetor$Lat.Long$Lat)+Range.Exp))
#'
#' RpraetorDataMat <- Array2Mat(Rpraetor$LMs)
#'
#' rThres <- IDbyDistanceRawDataCCV(LatLongs = Rpraetor$Lat.Long,
#'               RefData = RpraetorDataMat,
#'               ShapeData=TRUE,
#'               ShapeDim=2,
#'               DistMethod= "Proc",
#'               Verbose = TRUE,
#'               ProvConfidence = .95,
#'               PrintProg = FALSE,
#'               Method= 'Spearman')
#'
#' @export


IDbyDistanceRawDataCCV <- function(LatLongs, RefData, ShapeData=TRUE, ShapeDim=2, DistMethod=c("Euc", "Proc"), Verbose=TRUE, PrintProg=TRUE, ProvConfidence=0.95, Method=c('Spearman', 'Pearson')){

  UserInputAssessment(LatLongs, RefData, Method, RefDistMat = 'skip', DistVec = 'skip')

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
    ShapeDistMat <- ProcDistanceTable(Mat2Array(TotalShape, LMdim=ShapeDim))
  } else if (DistMethod=="Proc" && ShapeData==FALSE){
    stop("Error: Procrustes distance selected, but ShapeData argument is set to FALSE")
  }



  #creating an empty object to be populated by results
  CoordsHeat <- matrix(NA, nrow = dim(RefData)[1], ncol = 3)

  for (i in 1:dim(TotalShape)[1]){
    #i <- 1


    ShapeDist <- ShapeDistMat[i,-i]





    #This generates all the distances from the point on the map to all the specimen locations
    GeographicDist <- GeoDist2Point(RefLatLongs = LatLongs[-i,], TargetLatLong = LatLongs[i,])


    if (Method=='Spearman'){
      #running the correlation to generate r
      CorRes <- suppressWarnings(stats::cor.test(x = ShapeDist, y = GeographicDist, method = "spearman"))
    } else if (Method=='Pearson'){
      CorRes <- suppressWarnings(stats::cor.test(x = ShapeDist, y = GeographicDist, method = "pearson"))
    }


    #combining the individual results and organising them
    results <- c(LatLongs$Lats[i], LatLongs$Longs[i], CorRes$estimate)

    #adding the results of this round to the previous results
    CoordsHeat[i,] <- results

    if (PrintProg==TRUE){
      svMisc::progress(value = i, max.value = dim(TotalShape)[1], progress.bar = FALSE)
      Sys.sleep(0.01)

    }
  }

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

