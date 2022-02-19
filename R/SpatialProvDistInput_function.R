

#' Spatially provenance a specimen from distance data
#'
#' This function takes the distances from an unknown specimen to all reference specimens
#' and uses these distances to calculate a likely spatial provenance. Note that this
#' procedured can only be applied to one unknown specimen at a time. This distance input function
#' allows the user to input any dissimilarity or distance desired and as appropriate for the data.
#' The function has two applications either: calculating a specimens' provenance, or alternatively
#' it can be used to calculate the minimum correlation coefficient needed to correctly identify a
#' known specimen at its true collection location. The second application of this function can work
#' as a correct cross-validation process if looped, but see \code{\link{IDbyDistanceDistInputCCV}} function which
#' does this automatically in a leave-one-out process. However, if a corss-validation process that
#' removed more than one specimen from the reference dataset at a time is required then it is advised
#' that this be applied using this function.
#' @param DistDataVec is a vector of distances from each reference specimen to the specimen of interest (either a specimen with unknown provenance or another reference specimen that the user is interested in validating the provenance of)
#' @inheritParams IDbyDistanceRawData
#' @return If Verbose is TRUE then a dataframe of all values for every sampled grid reference is returned. If Verbose is FALSE then only those grid references with the highest correlation values are returned.
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
#' @examples
#' RatDistMat <- ProcDistanceTable(Rpraetor$LMs)
#' Range.Exp <- .5
#' Long.Range <- c(floor(min(Rpraetor$Lat.Long$Long))
#'               -Range.Exp,ceiling(max(Rpraetor$Lat.Long$Long)+Range.Exp))
#' Lat.Range <- c(floor(min(Rpraetor$Lat.Long$Lat))
#'               -Range.Exp,ceiling(max(Rpraetor$Lat.Long$Lat)+Range.Exp))
#' rThres <- IDbyDistanceDistInputCCV(LatLongs = Rpraetor$Lat.Long,
#'               DistDataMat = RatDistMat,
#'               Verbose = TRUE,
#'               ProvConfidence = .95,
#'               PrintProg = FALSE,
#'               Method = 'Spearman')
#'
#'
#' R.Samp <- c(12, 42)
#' IDbyDistanceDistInput(LatLongs = Rpraetor$Lat.Long[-1,],
#'                    DistDataVec = RatDistMat[-1,1],
#'                    LongRange = Long.Range,
#'                    LatRange = Lat.Range,
#'                    RangeSamp = R.Samp,
#'                    Verbose = FALSE,
#'                    Validate = FALSE,
#'                    PlotValCor = rThres$`Provenancing.Correlation.95%.Confidence`,
#'                    PlotProv = TRUE,
#'                    Method = 'Spearman')
#'
#'
#' @author Ardern Hulme-Beaman
#' @export


IDbyDistanceDistInput <- function(LatLongs, DistDataVec, LongRange, LatRange, RangeSamp=10, Verbose=TRUE, PrintProg=FALSE, Validate= FALSE, ValidLatLongs, PlotRes=TRUE, HeatHue= c(.15, 1), TileSize=2, PlotProv=FALSE, PlotValCor, Method=c('Pearson', 'Spearman'), PacificCent=FALSE){


  UserInputAssessment(LatLongs, RefData = NULL, RefDistMat = NULL, DistVec = DistDataVec, Method)


  #making LatLongs a dataframe
  LatLongs <- as.data.frame(LatLongs)
  colnames(LatLongs) <- c("Lats", "Longs")


  ShapeDist <- DistDataVec


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
    GeographicDist <- GeoDist2Point(RefLatLongs = LatLongs, TargetLatLong = ValidLatLongs)


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
    if (PacificCent==TRUE){
      MidRange <- seq(LongRange[1], LongRange[2]+360, by = LongRangeSteps*-1)
      #MidRange[which(MidRange>=180)] <- MidRange[which(MidRange>=180)]-360
      Longways <- c(LongRange[1]+LongRangeSteps, MidRange, LongRange[2]+(LongRangeSteps*-1)+360)
    } else {
      Longways <- c(LongRange[1]-LongRangeSteps, seq(LongRange[1], LongRange[2], by = LongRangeSteps), LongRange[2]+LongRangeSteps)
    }


    Latways <- c(LatRange[1]-LatRangeSteps, seq(LatRange[1], LatRange[2], by = LatRangeSteps), LatRange[2]+LatRangeSteps)

    if(sum(Latways<c(-90))>0){
      Latways <- Latways[-which(Latways<c(-90))]
    }
    if(sum(Latways>c(90))>0){
      Latways <- Latways[-which(Latways>c(90))]
    }



    #carrying out iterative analyses across the geographic region defined by LongRange and LatRange
    for (i in Longways){
      for (j in Latways){
        #i <- Longways[1]
        #j <- Latways[1]

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

  if (PacificCent==TRUE){
    PlottingMap <- "mapdata::world2Hires"
  } else {
    PlottingMap <- "world"
  }

  #if there is not a validating the data then we this means
  #we either have the validation result from a previous analyses and we can plot it
  #or we don't and we're not yet interested in it
  if (Validate==FALSE){
    if (PlotRes==TRUE){
      maps::map(PlottingMap, xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="black", bg=graphics::par(bg="white"))


      #creating colour scale from max and min slope based on which variable we're using
      CoordsHeatNum <- chr2nu(CoordsHeat$Cor)
      OriginLoc <- CoordsHeat[which(CoordsHeatNum==max(CoordsHeatNum)),]


      CoordsHeatscaled <- (CoordsHeatNum-min(CoordsHeatNum))/(max(CoordsHeatNum)-min(CoordsHeatNum))

      CoordsHeats <- grDevices::hsv(h = HeatHue[1], v = 1, s = CoordsHeatscaled, alpha = HeatHue[2])

      #plotting the correlations
      graphics::points(x = as.character(CoordsHeat$Longs), y = as.character(CoordsHeat$Lats), pch=15, col=CoordsHeats,  cex=TileSize)

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

      maps::map(PlottingMap, xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="black", bg=graphics::par(bg="white"), add=T)

      graphics::points(x = as.character(LatLongs$Longs), y = as.character(LatLongs$Lats), pch=23, bg='orange',  cex=1)
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



#' Spatial Provenancing Correct Cross-Validation calculation from distance data
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
#' @author Ardern Hulme-Beaman
#'
#' RatDistMat <- ProcDistanceTable(Rpraetor$LMs)
#'
#' Range.Exp <- .5
#'
#' Long.Range <- c(floor(min(Rpraetor$Lat.Long$Long))
#'               -Range.Exp,ceiling(max(Rpraetor$Lat.Long$Long)+Range.Exp))
#'
#' Lat.Range <- c(floor(min(Rpraetor$Lat.Long$Lat))
#'               -Range.Exp,ceiling(max(Rpraetor$Lat.Long$Lat)+Range.Exp))
#'
#' rThres <- IDbyDistanceDistInputCCV(LatLongs = Rpraetor$Lat.Long,
#'               DistDataMat = RatDistMat,
#'               Verbose = TRUE,
#'               ProvConfidence = .95,
#'               PrintProg = FALSE,
#'               Method = 'Spearman')
#'
#'
#' @export


IDbyDistanceDistInputCCV <- function(LatLongs, DistDataMat, Verbose=TRUE, PrintProg=TRUE, ProvConfidence=0.95, Method=c('Pearson', 'Spearman')){


  UserInputAssessment(LatLongs, RefData = NULL, DistVec = NULL, RefDistMat = DistDataMat, Method)

  #making LatLongs a dataframe
  LatLongs <- as.data.frame(LatLongs)
  colnames(LatLongs) <- c("Lats", "Longs")


  #creating an empty object to be populated by results
  CoordsHeat <- matrix(NA, nrow = dim(DistDataMat)[1], ncol = 3)

  for (i in 1:dim(DistDataMat)[1]){
    #i <- 1

    ShapeDist <- DistDataMat[-i,i]


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
      svMisc::progress(value = i, max.value = dim(DistDataMat)[1], progress.bar = FALSE)
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

