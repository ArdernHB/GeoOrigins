

#' Returns the geographic provenance of a specimen
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
#' @param Dist.data.vec is a vector of distances from each reference specimen to the specimen of interest (either a specimen with unknown provenance or another reference specimen that the user is interested in validating the provenance of)
#' @inheritParams IDbyDistance.RawData
#' @return If verbose is TRUE then a dataframe of all values for every sampled grid reference is returned. If verbose is FALSE then only those grid references with the highest correlation values are returned.
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
#' @keywords Spatial provenancing
#' @keywords Spatial identification
#' @author Ardern Hulme-Beaman
#' @export


IDbyDistance.DistInput <- function(Lat.Longs=data.frame(Lats=c(), Longs=c()), Dist.data.vec=Dist, LongRange=c(0,0), LatRange=c(0,0), Range.Samp=10, verbose=TRUE, print.prog=FALSE, Validate= FALSE, Valid.LatLongs=c(Lat, Long), Plot.Res=TRUE, HeatHue= c(.1, 1), Tile.Size=2, plot.Prov=FALSE, plot.Val.cor=as.numeric(X)){

  #making Lat.Longs a dataframe
  Lat.Longs <- as.data.frame(Lat.Longs)
  colnames(Lat.Longs) <- c("Lats", "Longs")


  Shape.Dist <- Dist.data.vec


  #creating an empty object to be populated by results
  CoordsHeat <- NULL

  #this function has two operations
  #either it can come up with the correlation value at a specific point (Validate=TRUE)
  #or it can come up with correlation values across the entire map
  #here we have the first option which is useful for the validation process
  #by doing this first option for all the specimens (in a loop) we can build a distribution
  #of the correlation values that will correctly cover the specimens true location
  #but see IDbyDistance.RawData.ccv function for looping to be done automatically
  if (Validate==TRUE){

    Valid.LatLongs <- as.data.frame(Valid.LatLongs)
    colnames(Valid.LatLongs) <- c("Lats", "Longs")

    #This generates all the distances from the point on the map to all the specimen locations
    Geographic.Dists <- Geo.Dist2Point(RefLatLongs = Lat.Longs, TargetLatLong = c(Valid.LatLongs$Lats, Valid.LatLongs$Longs))


    #running the correlation to generate r
    CorRes <- stats::cor.test(x = Shape.Dist, y = Geographic.Dists)


    #combining the individual results and organising them
    results <- c(Valid.LatLongs$Lats, Valid.LatLongs$Longs, CorRes$estimate)

    #adding the results of this round to the previous results
    CoordsHeat <- rbind(CoordsHeat, results)
  } else if (Validate==FALSE){
    #this is the second of the above options and not for validation
    #here the function goes through the entire map and generated correlation values
    #the level of detail is set by "Range.Samp"
    #for example if Range.Samp is set to 1 the following procedure
    #will generate correlation values at every degree of the map
    #if the Range.Samp is set to 2 it will generate for every other degree
    #this can speed things up for quick and dirty tests

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


    #carrying out iterative analyses across the geographic region defined by LongRange and LatRange
    for (i in Longways){
      for (j in Latways){
        coord <- c(j, i)

        Geographic.Dists <- Geo.Dist2Point(RefLatLongs = Lat.Longs, TargetLatLong = coord)

        #running the correlation to generate r
        CorRes <- cor.test(x = Shape.Dist, y = Geographic.Dists)

        #combining the individual results and organising them
        results <- c(as.character(j),i, CorRes$estimate)

        #adding the results of this round to the previous results
        CoordsHeat <- rbind(CoordsHeat, results)
      }

      if (print.prog==TRUE){
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
    if (Plot.Res==TRUE){
      maps::map("world", xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="black", bg=par(bg="white"))


      #creating colour scale from max and min slope based on which variable we're using
      CoordsHeatNum <- chr2nu(CoordsHeat$Cor.Slope)
      OriginLoc <- CoordsHeat[which(CoordsHeatNum==max(CoordsHeatNum)),]


      CoordsHeatscaled <- (CoordsHeatNum-min(CoordsHeatNum))/(max(CoordsHeatNum)-min(CoordsHeatNum))

      CoordsHeats <- hsv(h = HeatHue[1], v = 1, s = CoordsHeatscaled, alpha = HeatHue[2])

      #plotting the correlations
      points(x = as.character(CoordsHeat$Long), y = as.character(CoordsHeat$Lat), pch=15, col=CoordsHeats,  cex=Tile.Size)

      #here if we want to plot a polygon of the region that
      #approximates the region the specimens came from (with whatever level of confidence we have selected)
      if (plot.Prov==TRUE){
        Select.95.conf <- which(chr2nu(CoordsHeat$Cor.Slope)>plot.Val.cor)

        ApproxOrigin <- CoordsHeat[Select.95.conf,]

        contour.95 <-  Construct_contour(ApproxOrigin)
        polycol <- hsv(h = HeatHue[1], s = 1, v = .8, alpha = HeatHue[2])

        polygon(contour.95$x, contour.95$y, col=NA, border=polycol, lwd=2)

      }

      maps::map("world", xlim=c(min(Longways), max(Longways)), ylim=c(min(Latways), max(Latways)), interior=FALSE, col="black", bg=par(bg="white"), add=T)

      points(x = as.character(Lat.Longs$Long), y = as.character(Lat.Longs$Lat), pch=23, bg='orange',  cex=1)
      maps::map.axes()
    }
  }


  if (verbose==TRUE|Validate==TRUE){
    return(CoordsHeat)
  } else {
    OriginLocCor <- CoordsHeat[which(CoordsHeat$Cor.Slope==max(chr2nu(CoordsHeat$Cor.Slope), na.rm = TRUE)),]

    return(list(Cor=OriginLocCor))
  }


}



#' Spatial Provenancing Correct Cross-Validation calculation from distance data
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
#'  @details This method also makes use of the \code{cor.test} function from the \code{stats} package. When the \code{print.prog} is set to TRUE, the \code{progress} function of the \code{svMisc} package is used.
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


IDbyDistance.DistInput.CCV <- function(Lat.Longs=data.frame(Lats=c(), Longs=c()), Dist.data.mat=DistMat, verbose=TRUE, print.prog=TRUE, Prov.Confidence=0.95){

  #making Lat.Longs a dataframe
  Lat.Longs <- as.data.frame(Lat.Longs)
  colnames(Lat.Longs) <- c("Lats", "Longs")


  #creating an empty object to be populated by results
  CoordsHeat <- NULL

  for (i in 1:dim(total.shape)[1]){
    #i <- 1

    Shape.Dist <- Dist.data.mat[-i,i]


    #This generates all the distances from the point on the map to all the specimen locations
    Geographic.Dists <- Geo.Dist2Point(RefLatLongs = Lat.Longs[-i,], TargetLatLong = Lat.Longs[i,])


    #running the correlation to generate r
    CorRes <- stats::cor.test(x = Shape.Dist, y = Geographic.Dists)


    #combining the individual results and organising them
    results <- c(Valid.LatLongs$Lats, Valid.LatLongs$Longs, CorRes$estimate)

    #adding the results of this round to the previous results
    CoordsHeat <- rbind(CoordsHeat, results)

    if (print.prog==TRUE){
      svMisc::progress(value = i, max.value = dim(total.shape)[1], progress.bar = FALSE)
      Sys.sleep(0.01)

    }
  }

  #converting the results into a data frame
  CoordsHeat <- as.data.frame(CoordsHeat, row.names = 1:dim(CoordsHeat)[1])

  #naming the variables
  names(CoordsHeat) <- c("Lats", "Longs", "Cor")

  ProvCor <- quantile(CoordsHeat$Cor, 1-Prov.Confidence)

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

